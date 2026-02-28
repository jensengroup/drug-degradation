import csv
import json
from collections import defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor, rdMolDescriptors
from tqdm import tqdm

from rxnnet.thermo import Thermochemistry
from rxnnet.config import Config

rdDepictor.SetPreferCoordGen(True)

# Constants
HARTREE_TO_KCALMOL = 627.5094740631


def hartree2kcalmol(energy: float) -> float:
    """Convert energy from Hartree to kcal/mol."""
    return energy * HARTREE_TO_KCALMOL


class QMResult:
    def __init__(self, l1, l2):
        self.l1 = l1
        self.l2 = l2
        self.validate()
        self.thermo_calculator = Thermochemistry(
            *self._get_ac(self.l1),
            self.vibs,
        )

    @classmethod
    def from_files(cls, file1, file2):
        with open(file1) as f1:
            data1 = json.load(f1)
        with open(file2) as f2:
            data2 = json.load(f2)
        return cls(data1, data2)

    def validate(self):
        try:
            assert "Hessian" in self.l1["Geometries"][-1], "L1 file missing Hessian"
            assert self._coords_match, "L1 and L2 coordinates do not match"
            assert "Single_Point_Data" in self.l2["Geometries"][-1], (
                "L2 file missing Single_Point_Data"
            )
            assert self.l2["Geometries"][-1]["Single_Point_Data"]["Converged"], (
                "L2 file Single_Point_Data not converged"
            )
        except Exception as e:
            raise ValueError(e)

    @staticmethod
    def _get_ac(data):
        tmp = data["Geometries"][-1]["Geometry"]["Coordinates"]["Cartesians"]
        atoms = np.array([a[0] for a in tmp])
        coords = np.array([c[1:] for c in tmp], dtype=float) * 0.529177249
        return atoms, coords

    @property
    def _coords_match(self):
        coords_l1 = self._get_ac(self.l1)[1]
        coords_l2 = self._get_ac(self.l2)[1]
        return np.allclose(coords_l1, coords_l2)

    @property
    def vibs(self):
        vibs = self.l1["Geometries"][-1]["THERMOCHEMISTRY_Energies"][-1]["FREQ"]
        vibs = np.array(
            [
                f[0]
                * self.l1["Geometries"][-1]["THERMOCHEMISTRY_Energies"][-1][
                    "freqScalingFactor"
                ]
                for f in vibs
            ]
        )
        return vibs[vibs > 0]

    def get_energy(self, T=298.15, p_atm=1.0, c_M=None):
        # calculate gibbs correction
        gibbs_corr = self.thermo_calculator.get_contributions(
            T=T, p_atm=p_atm, c_M=c_M
        )["gibbs_correction_Eh"]
        # electronic energy from L2
        e_elec = self.l2["Geometries"][-1]["Single_Point_Data"]["FinalEnergy"]
        return float(e_elec + gibbs_corr)


@dataclass
class MoleculeNode:
    """A node in the reaction network representing a molecular species."""

    idx: int
    smiles: str = ""
    mol: Optional[Chem.Mol] = None
    energy: float = 0.0  # Electronic energy in Hartree
    charge: int = 0
    weight: float = 0.0
    svg: str = ""
    origin_type: str = ""
    origin_ids: List[int] = field(default_factory=list)
    sibling_ids: List[int] = field(default_factory=list)


@dataclass
class Reaction:
    reactant_ids: List[int]
    product_ids: List[int]  # First is primary, rest are eliminated
    rxn_type: str
    count: int = 1
    reaction_smiles: str = ""
    ts_mol: Optional[Chem.Mol] = None
    ts_energy: Optional[float] = None
    ts_molblock: str = ""

    @property
    def id(self) -> str:
        r_str = "_".join(map(str, sorted(self.reactant_ids)))
        p_str = "_".join(map(str, sorted(self.product_ids)))
        return f"{r_str}->{p_str}:{self.rxn_type}"

    @property
    def primary_reactant(self) -> int:
        return self.reactant_ids[0] if self.reactant_ids else -1

    @property
    def primary_product(self) -> int:
        return self.product_ids[0] if self.product_ids else -1

    @property
    def eliminated_products(self) -> List[int]:
        return self.product_ids[1:] if len(self.product_ids) > 1 else []

    @property
    def is_multi_product(self) -> bool:
        return len(self.product_ids) > 1


class ReactionNetwork:
    # Constants
    R = 1.98720425864083 / 1000  # Gas constant in kcal/(mol*K)
    PKA_SLOPE = 0.43496246
    PKA_INTERCEPT = -114.66798804

    STANDARD_STATES = {"[H]O[H]": 55.34}

    def __init__(
        self,
        network_dir: str | Path,
        substrate_id: int = 1,
        pH: float = 7.0,
        temperature: float = 313.15,
    ):
        self.config = Config(network_dir)
        self.substrate_id = substrate_id
        self.pH = pH
        self.temperature = temperature

        self.nodes: Dict[int, MoleculeNode] = {}
        self.reactions: List[Reaction] = []

        self._reactant_index: Dict[int, List[Reaction]] = {}
        self._product_index: Dict[int, List[Reaction]] = {}
        self._edge_map: Dict[str, Reaction] = {}

        # Pathway data
        self.pathways: Dict[int, List[Dict]] = {}

        # Load
        self._load_network()

    def _load_network(self) -> None:
        print(f"Loading network from {self.config.network_dir}")
        self._load_nodes()
        self._load_reactions()
        self._build_indexes()
        print(f"Loaded {len(self.nodes)} nodes and {len(self.reactions)} reactions")

    def _load_nodes(self) -> None:
        nodes_csv = self.config.network_dir / "nodes.csv"
        if not nodes_csv.exists():
            return

        smiles_map = {}
        with open(nodes_csv, "r") as f:
            for row in csv.DictReader(f):
                smiles_map[int(row["idx"])] = row.get("canonical_smiles", "")

        node_files = list(self.config.node_data.glob("*.json"))
        for node_file in tqdm(node_files, desc="Loading nodes"):
            try:
                with open(node_file, "r") as f:
                    data = json.load(f)

                idx = data["idx"]
                mol_block = data.get("data", {}).get("gfn2-xtb", "")
                mol = (
                    Chem.MolFromMolBlock(mol_block, removeHs=False)
                    if mol_block
                    else None
                )
                props = self._parse_sdf_properties(mol_block) if mol_block else {}

                origin_ids = []
                sibling_ids = []
                if "origin-ids" in props:
                    try:
                        origin_ids = eval(props["origin-ids"])
                    except KeyError:
                        pass
                if "sibling_labels" in props:
                    try:
                        sibling_ids = eval(props["sibling_labels"])
                    except KeyError:
                        pass

                # get energy of node
                level = "quick"
                f_l1 = self.config.qm_data / f"{idx}-freq-{level}.json"
                f_l2 = self.config.qm_data / f"{idx}-sp-{level}.json"
                if not f_l1.is_file() or not f_l2.is_file():
                    energy = 0.0
                else:
                    standard_state = self.STANDARD_STATES.get(
                        smiles_map.get(idx, ""), 1.0
                    )
                    result = QMResult.from_files(f_l1, f_l2)
                    energy = result.get_energy(
                        T=self.temperature, c_M=standard_state, p_atm=None
                    )

                self.nodes[idx] = MoleculeNode(
                    idx=idx,
                    smiles=smiles_map.get(idx, ""),
                    mol=mol,
                    energy=energy,
                    charge=Chem.rdmolops.GetFormalCharge(mol) if mol else 0,
                    weight=round(rdMolDescriptors.CalcExactMolWt(mol)) if mol else 0,
                    svg=self._generate_svg(mol) if mol else "",
                    origin_type=props.get("origin-type", ""),
                    origin_ids=origin_ids,
                    sibling_ids=sibling_ids,
                )
            except Exception as e:
                print(f"Error loading {node_file}: {e}")

    def _load_reactions(self) -> None:
        reactions_csv = self.config.network_dir / "reactions.csv"
        if not reactions_csv.exists():
            return

        with open(reactions_csv, "r") as f:
            for row in csv.DictReader(f):
                try:
                    reactant_ids = json.loads(row["reactant_ids"])
                    product_ids = json.loads(row["product_ids"])

                    # Sort products by size (largest first)
                    if len(product_ids) > 1:
                        product_ids = sorted(
                            product_ids,
                            key=lambda x: (
                                self.nodes[x].mol.GetNumAtoms()
                                if x in self.nodes and self.nodes[x].mol
                                else 0
                            ),
                            reverse=True,
                        )

                    rxn = Reaction(
                        reactant_ids=reactant_ids,
                        product_ids=product_ids,
                        rxn_type=row.get("rxn_type", "unknown"),
                        count=int(row.get("count", 1)),
                        reaction_smiles=row.get("reaction_smiles", ""),
                    )

                    # Load TS if available
                    if len(reactant_ids) == 1 and product_ids:
                        ts_file = (
                            self.config.qm_data
                            / f"ts-{reactant_ids[0]}-{product_ids[0]}.sdf"
                        )
                        if ts_file.exists():
                            suppl = Chem.SDMolSupplier(
                                str(ts_file), removeHs=False, sanitize=False
                            )
                            ts_mol = next(iter(suppl), None) if suppl else None
                            if ts_mol:
                                rxn.ts_mol = ts_mol
                                rxn.ts_molblock = Chem.MolToMolBlock(ts_mol)
                                if ts_mol.HasProp("electronic_energy"):
                                    rxn.ts_energy = float(
                                        ts_mol.GetProp("electronic_energy")
                                    )

                    self.reactions.append(rxn)
                except Exception as e:
                    print(f"Error loading reaction: {e}")

    def _build_indexes(self) -> None:
        for rxn in self.reactions:
            for r_id in rxn.reactant_ids:
                self._reactant_index.setdefault(r_id, []).append(rxn)
            for p_id in rxn.product_ids:
                self._product_index.setdefault(p_id, []).append(rxn)
            if rxn.reactant_ids and rxn.product_ids:
                self._edge_map[f"{rxn.primary_reactant}-{rxn.primary_product}"] = rxn

    @staticmethod
    def _parse_sdf_properties(mol_block: str) -> Dict[str, str]:
        props = {}
        lines = mol_block.split("\n")
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith(">") and "<" in line:
                start, end = line.index("<") + 1, line.index(">", line.index("<"))
                prop_name = line[start:end]
                i += 1
                values = []
                while (
                    i < len(lines)
                    and lines[i].strip()
                    and not lines[i].strip().startswith(">")
                ):
                    values.append(lines[i].strip())
                    i += 1
                props[prop_name] = "\n".join(values)
            else:
                i += 1
        return props

    def _generate_svg(self, mol: Chem.Mol) -> str:
        if not mol:
            return ""

        draw_mol = Chem.RemoveHs(mol)
        for a in draw_mol.GetAtoms():
            a.SetAtomMapNum(0)
        draw_mol.RemoveAllConformers()
        return Draw.MolsToGridImage(
            [draw_mol], molsPerRow=1, subImgSize=(300, 250), useSVG=True
        )

    def get_node(self, idx: int) -> Optional[MoleculeNode]:
        return self.nodes.get(idx)

    def get_reactions_from(self, node_idx: int) -> List[Reaction]:
        return self._reactant_index.get(node_idx, [])

    def get_reaction_between(self, r_id: int, p_id: int) -> Optional[Reaction]:
        return self._edge_map.get(f"{r_id}-{p_id}")

    def get_products(self, node_idx: int) -> Set[int]:
        return {rxn.primary_product for rxn in self.get_reactions_from(node_idx)}

    def compute_pathways(self, max_extra_steps: int = 3, max_paths: int = 100) -> None:
        """Compute all pathways from substrate to all reachable nodes."""
        if self.substrate_id is None:
            print("Warning: No substrate_id set")
            return

        # Build adjacency
        graph = defaultdict(list)
        for rxn in self.reactions:
            for r_id in rxn.reactant_ids:
                graph[r_id].append(rxn.primary_product)

        # BFS for shortest distances
        distances = {self.substrate_id: 0}
        queue = deque([self.substrate_id])
        while queue:
            curr = queue.popleft()
            for neighbor in graph[curr]:
                if neighbor not in distances:
                    distances[neighbor] = distances[curr] + 1
                    queue.append(neighbor)

        # DFS for all paths
        all_paths = defaultdict(list)

        def dfs(node, path, dist):
            for target, target_dist in distances.items():
                if target != self.substrate_id:
                    max_dist = target_dist + max_extra_steps
                    if (
                        node == target
                        and dist <= max_dist
                        and len(all_paths[target]) < max_paths
                    ):
                        all_paths[target].append(path.copy())

            for neighbor in graph[node]:
                if neighbor not in path:
                    path.append(neighbor)
                    dfs(neighbor, path, dist + 1)
                    path.pop()

        dfs(self.substrate_id, [self.substrate_id], 0)

        # Compute profiles
        for target, paths in tqdm(all_paths.items(), desc="Computing profiles"):
            self.pathways[target] = []
            for path in paths:
                profile = self._compute_profile(path)
                barrier = self._max_barrier(profile)
                self.pathways[target].append(
                    {
                        "path": path,
                        "profile": profile,
                        "barrier_height": barrier,
                        "path_length": len(path) - 1,
                    }
                )

        reachable = sum(1 for p in self.pathways.values() if p)
        print(f"Computed pathways for {reachable} reachable nodes")

    def _compute_profile(self, path: List[int]) -> List[Dict]:
        """Compute energy profile for a path with pH corrections."""
        if len(path) < 2:
            return []

        T = self.temperature
        eliminated_offset = 0.0
        proton_offset = 0.0
        proton_type = ""
        profile = []
        coord = 0.0

        first_energy = hartree2kcalmol(self.nodes[path[0]].energy)
        profile.append(
            {
                "node_id": path[0],
                "energy": first_energy,
                "type": "minimum",
                "reaction_coordinate": coord,
                "step": 0,
            }
        )
        prev_energy = first_energy

        for i in range(1, len(path)):
            rxn = self.get_reaction_between(path[i - 1], path[i])
            if not rxn:
                continue

            node = self.nodes[path[i]]
            prev_node = self.nodes[path[i - 1]]
            edge_type = rxn.rxn_type.lower()

            # Transition state
            if rxn.ts_mol and rxn.ts_energy:
                ts_energy = prev_energy + hartree2kcalmol(
                    rxn.ts_energy - prev_node.energy
                )
                coord += 0.5
                profile.append(
                    {
                        "node_id": f"TS{i}",
                        "energy": ts_energy,
                        "type": "transition_state",
                        "reaction_coordinate": coord,
                        "step": i - 0.5,
                    }
                )
                coord += 0.5
            else:
                coord += 1.0

            # Eliminated products
            for elim_id in rxn.eliminated_products:
                if elim_id in self.nodes:
                    eliminated_offset += hartree2kcalmol(self.nodes[elim_id].energy)

            # Protonation/deprotonation corrections
            if edge_type == "protonation" or edge_type == "protonate":
                if proton_type == "deprotonation":
                    proton_type, proton_offset = "", 0.0
                else:
                    dG = hartree2kcalmol(node.energy) - prev_energy
                    pka = self.PKA_SLOPE * (-dG) + self.PKA_INTERCEPT
                    proton_offset = (
                        prev_energy
                        - hartree2kcalmol(node.energy)
                        + self.R * T * np.log(10) * (self.pH - pka)
                    )
                    proton_type = "protonation"
            elif edge_type == "deprotonation" or edge_type == "deprotonate":
                if proton_type == "protonation":
                    proton_type, proton_offset = "", 0.0
                else:
                    dG = hartree2kcalmol(node.energy) - prev_energy
                    pka = self.PKA_SLOPE * dG + self.PKA_INTERCEPT
                    proton_offset = (
                        prev_energy
                        - hartree2kcalmol(node.energy)
                        + self.R * T * np.log(10) * (pka - self.pH)
                    )
                    proton_type = "deprotonation"

            node_energy = (
                hartree2kcalmol(node.energy) + eliminated_offset + proton_offset
            )
            prev_energy = node_energy

            profile.append(
                {
                    "node_id": path[i],
                    "energy": node_energy,
                    "type": "minimum",
                    "reaction_coordinate": coord,
                    "step": i,
                    "edge_type": edge_type,
                }
            )

        # Normalize to substrate = 0
        base = profile[0]["energy"]
        for p in profile:
            p["energy"] -= base

        return profile

    def _max_barrier(self, profile: List[Dict]) -> float:
        if not profile:
            return 0.0
        energies = [p["energy"] for p in profile]
        rolling_min = float("inf")
        max_barrier = 0.0
        for e in energies:
            rolling_min = min(rolling_min, e)
            max_barrier = max(max_barrier, e - rolling_min)
        return max_barrier

    def get_best_pathway(self, node_id: int) -> Optional[Dict]:
        """Get lowest-barrier pathway to a node."""
        if node_id not in self.pathways or not self.pathways[node_id]:
            return None
        return min(self.pathways[node_id], key=lambda p: p["barrier_height"])

    def to_visualization_data(self) -> Dict[str, Any]:
        """Export data for JavaScript visualization."""
        # Node maps
        svg_map = {idx: n.svg for idx, n in self.nodes.items()}
        energy_map = {idx: n.energy for idx, n in self.nodes.items()}
        charge_map = {idx: n.charge for idx, n in self.nodes.items()}
        weight_map = {idx: n.weight for idx, n in self.nodes.items()}

        # Relative energies
        rel_energy_map = {}
        if self.substrate_id and self.substrate_id in self.nodes:
            sub_e = self.nodes[self.substrate_id].energy
            rel_energy_map = {
                idx: hartree2kcalmol(n.energy - sub_e) for idx, n in self.nodes.items()
            }

        # Edge data
        edge_data = []
        for rxn in self.reactions:
            e = {
                "begin": rxn.primary_reactant,
                "end": rxn.primary_product,
                "type": rxn.rxn_type,
                "count": rxn.count,
                "smaller_products": rxn.eliminated_products,
            }
            if rxn.ts_mol:
                e["ts"] = rxn.ts_molblock
                e["ts_energy"] = rxn.ts_energy
            edge_data.append(e)

        # Best paths
        best_energy = {}
        best_barrier = {}
        best_path = {}
        for node_id, pathways in self.pathways.items():
            if pathways:
                best = min(pathways, key=lambda p: p["barrier_height"])
                best_energy[node_id] = (
                    best["profile"][-1]["energy"] if best["profile"] else 0
                )
                best_barrier[node_id] = best["barrier_height"]
                best_path[node_id] = best["path"]

        # Stereoisomer groups
        stereo_groups = {}
        for idx, node in self.nodes.items():
            if node.sibling_ids:
                group = [idx] + node.sibling_ids
                key = min(group)
                if key not in stereo_groups:
                    stereo_groups[key] = group

        return {
            "svgMap": svg_map,
            "molEnergyMap": energy_map,
            "chargeMap": charge_map,
            "weightMap": weight_map,
            "relEnergyMap": rel_energy_map,
            "originalEdgeData": edge_data,
            "substrateId": self.substrate_id,
            "temperature": self.temperature,
            "pH": self.pH,
            "bestPathEnergyMap": best_energy,
            "bestPathBarrierMap": best_barrier,
            "bestPathMap": best_path,
            "stereoisomerGroups": stereo_groups,
            "pathways": {str(k): v for k, v in self.pathways.items()},
        }

    def __repr__(self) -> str:
        return (
            f"ReactionNetwork(nodes={len(self.nodes)}, reactions={len(self.reactions)})"
        )
