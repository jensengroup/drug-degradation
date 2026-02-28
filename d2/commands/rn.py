"""
Network Visualization Module - Redesigned for simplicity and maintainability.

This module provides tools for visualizing chemical reaction networks with
a focus on clarity and maintainability. All calculations are performed
on the Python side, with JavaScript used only for visualization and interactivity.
"""

import shelve
from typing import Any, Dict, List, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from tooltoad.chemutils import hartree2kcalmol
from tqdm import tqdm

from d2.config import NetworkConfig
from d2.network import Node
from d2.utils import Reaction

rdDepictor.SetPreferCoordGen(True)


class ReactionNetwork:
    """Simplified reaction network class for loading data and analyzing
    pathways."""

    # Gas constant in kcal/(mol*K)
    R = 1.98720425864083 / 1000

    STANDARD_STATES = {"[H]O[H]": 55.34}

    def __init__(
        self,
        config: NetworkConfig,
        substrate_id: int,
        prop_name: str,
        pH: float = 7.0,
        level: str = "quick",
    ):
        """Initialize the reaction network and precompute all pathways.

        Args:
            config: The network configuration
            substrate_id: The ID of the substrate molecule
            prop_name: The property name to use for energy calculations
            pH: The pH value to use for energy calculations
        """
        self.config = config
        self.substrate_id = substrate_id
        self.prop_name = prop_name
        self.pH = pH
        self.energy_type = "G" if "gibbs" in prop_name else "E"

        self.nodes: Dict[int, Dict] = {}
        self.edges: List[Dict] = []
        self.node_labels: Dict[int, str] = {}
        self.edge_map: Dict[str, Dict] = {}

        # precomputed pathway data
        self.pathways: Dict[int, List[Dict]] = {} 

        self._profile_cache: Dict[tuple, List[Dict]] = {}

        with shelve.open(str(self.config.label_db)) as db:
            self.label2id = db.get("label2id", {})

        print("Loading reaction network data...")
        self._load_nodes(level=level)
        self._process_reactions()  
        self._calculate_relative_energies()

        print("Precomputing all pathways and reaction profiles...")
        self._precompute_all_pathways()

        print(f"Network loaded: {len(self.nodes)} nodes, {len(self.edges)} edges")
        reachable = sum(1 for pathways in self.pathways.values() if pathways)
        print(f"Precomputed pathways for {reachable} reachable nodes")

        ts_count = sum(1 for edge in self.edges if "ts_energy" in edge)
        print(f"Found transition states for {ts_count} reactions")

    def _precompute_all_pathways(
        self, max_extra_steps: int = 3, max_paths_per_node: int = 100
    ):
        """Two-phase approach: find paths quickly, then calculate profiles selectively."""
        from collections import defaultdict

        graph = defaultdict(list)

        for edge in self.edges:
            begin_id = edge["begin"]
            end_id = edge["end"]
            graph[begin_id].append(end_id)

        all_node_paths = self._find_all_paths_bulk(
            graph, max_extra_steps, max_paths_per_node
        )

        for target_node, paths in tqdm(all_node_paths.items(), desc="Processing nodes"):
            if not paths:
                self.pathways[target_node] = []
                continue

            # Calculate full profiles for all paths
            pathway_data = []
            for path in paths:
                try:
                    profile = self._calculate_reaction_profile(path, self.edge_map)
                    barrier = self._find_max_barrier(profile)
                    # see if there are smaller products somewhere and add that to the path data

                    pathway_data.append(
                        {
                            "path": path,
                            "profile": profile,
                            "barrier_height": barrier,
                            "path_length": len(path) - 1,
                        }
                    )
                except Exception as e:
                    print(f"Error calculating profile for path {path}: {e}")

            self.pathways[target_node] = pathway_data

    def _find_all_paths_bulk(self, graph, max_extra_steps, max_paths_per_node):
        """Find paths to all nodes in a single traversal - much more efficient."""
        from collections import defaultdict, deque

        # find shortest distances to all nodes
        shortest_distances = {self.substrate_id: 0}
        queue = deque([self.substrate_id])

        while queue:
            current = queue.popleft()
            current_dist = shortest_distances[current]

            for neighbor in graph[current]:
                if neighbor not in shortest_distances:
                    shortest_distances[neighbor] = current_dist + 1
                    queue.append(neighbor)

        # collect paths within tolerance using optimized DFS
        all_paths = defaultdict(list)

        def dfs_all_targets(current_node, path, current_dist):
            min_remaining_dist = float("inf")

            for target in shortest_distances:
                if target != self.substrate_id:
                    target_dist = shortest_distances[target] + max_extra_steps
                    remaining = target_dist - current_dist

                    if remaining >= 0:
                        min_remaining_dist = min(min_remaining_dist, remaining)

                        # Check if current node is a valid target
                        if (
                            current_node == target
                            and current_dist <= target_dist
                            and len(all_paths[target]) < max_paths_per_node
                        ):
                            all_paths[target].append(path.copy())

            if min_remaining_dist > 0:
                for neighbor in graph[current_node]:
                    if neighbor not in path:  # Avoid cycles
                        if neighbor in shortest_distances:
                            neighbor_min_dist = shortest_distances[neighbor]
                            if current_dist + 1 + neighbor_min_dist <= max(
                                shortest_distances[t] + max_extra_steps
                                for t in shortest_distances
                                if t != self.substrate_id
                            ):
                                path.append(neighbor)
                                dfs_all_targets(neighbor, path, current_dist + 1)
                                path.pop()

        dfs_all_targets(self.substrate_id, [self.substrate_id], 0)

        total_paths = sum(len(paths) for paths in all_paths.values())
        print(f"Found {total_paths} total paths to {len(all_paths)} nodes")

        return dict(all_paths)

    def _calculate_reaction_profile(
        self, path: List[int], edge_map: Dict[str, Dict]
    ) -> List[Dict]:
        """Calculate energy profile for a reaction path including transition
        states."""
        if len(path) < 2:
            return []

        # Constants for pKa calculations
        T = 313.15
        PKA_SLOPE = 0.43496246
        PKA_INTERCEPT = -114.66798804

        eliminated_offset = 0.0
        proton_offset = 0.0
        proton_type = ""
        eliminated = []
        profile = []
        energies = [hartree2kcalmol(self.nodes[path[0]]["energy"])]
        reaction_coordinate = 0.0

        
        profile.append(
            {
                "node_id": path[0],
                "energy": energies[-1],
                "type": "minimum",
                "reaction_coordinate": reaction_coordinate,
                "step": 0,
                "eliminated": eliminated.copy(),
            }
        )

        for i in range(1, len(path)):
            edge_key = f"{path[i-1]}-{path[i]}"
            if edge_key not in edge_map:
                continue

            edge = edge_map[edge_key]
            node = self.nodes[path[i]]
            edge_type = edge["type"].lower()
            # Add transition state if it exists
            if "ts" in edge and edge["ts"]:
                try:
                    ts_mol = Chem.MolFromMolBlock(edge["ts"])
                    ts_energy_raw = edge.get("ts_energy")
                    if ts_mol and ts_energy_raw:
                        # take the previous nodes processed energy and add the difference of the raw Ts and previous nodes raw energy
                        ts_energy = energies[-1] + hartree2kcalmol(
                            ts_energy_raw - self.nodes[path[i - 1]]["energy"]
                        )
                        reaction_coordinate += 0.5
                        profile.append(
                            {
                                "node_id": f"TS{i}",
                                "energy": ts_energy,
                                "type": "transition_state",
                                "reaction_coordinate": reaction_coordinate,
                                "step": i - 0.5,
                                "ts_data": edge["ts"],
                                "eliminated": eliminated.copy(),
                            }
                        )

                        reaction_coordinate += 0.5
                    else:
                        reaction_coordinate += 1.0
                except Exception as e:
                    print(f"Warning: Could not process TS for edge {edge_key}: {e}")
                    reaction_coordinate += 1.0
            else:
                reaction_coordinate += 1.0

            # Handle eliminated products
            if "smaller_products" in edge:
                for idx in edge["smaller_products"]:
                    eliminated_offset += hartree2kcalmol(self.nodes[idx]["energy"])
                    eliminated.append(idx)
            # Handle protonation/deprotonation with offset
            if edge_type == "protonation":
                if proton_type == "deprotonation":
                    # Neutralization
                    proton_type = ""
                    proton_offset = 0.0
                else:
                    dG = hartree2kcalmol(node["energy"]) - energies[-1]
                    pka = PKA_SLOPE * (-dG) + PKA_INTERCEPT
                    dG_correction = self.R * T * np.log(10) * (self.pH - pka)
                    proton_offset = (
                        energies[-1] - hartree2kcalmol(node["energy"]) + dG_correction
                    )
                    proton_type = node.get("origin-type", "").lower()

            elif edge_type == "deprotonation":
                if proton_type == "protonation":
                    # Neutralization
                    proton_type = ""
                    proton_offset = 0.0
                else:
                    dG = hartree2kcalmol(node["energy"]) - energies[-1]
                    pka = PKA_SLOPE * dG + PKA_INTERCEPT
                    dG_correction = self.R * T * np.log(10) * (pka - self.pH)
                    proton_offset = (
                        energies[-1] - hartree2kcalmol(node["energy"]) + dG_correction
                    )
                    proton_type = node.get("origin-type", "").lower()

            node_energy = (
                hartree2kcalmol(node["energy"]) + eliminated_offset + proton_offset
            )
            energies.append(node_energy)

            profile.append(
                {
                    "node_id": path[i],
                    "energy": node_energy,
                    "type": "minimum",
                    "reaction_coordinate": reaction_coordinate,
                    "step": i,
                    "edge_type": edge_type,
                    "eliminated": eliminated.copy(),
                }
            )

        energies = np.array([p["energy"] for p in profile])
        substrate_energy = energies[0]

        for i, point in enumerate(profile):
            point["energy"] = energies[i] - substrate_energy

        return profile

    def _find_max_barrier(self, profile: List[Dict]) -> float:
        """Find maximum energy barrier in a reaction profile."""
        if not profile:
            return 0.0

        energies = [p["energy"] for p in profile]
        rolling_min = float("inf")
        max_barrier = 0.0

        for energy in energies:
            if energy < rolling_min:
                rolling_min = energy
            barrier = energy - rolling_min
            max_barrier = max(max_barrier, barrier)

        return max_barrier

    def filter_pathways(
        self, min_count: Optional[int] = None, max_barrier: Optional[float] = None
    ) -> Dict[int, List[Dict]]:
        """Enhanced filtering with early termination and better sorting.

        Args:
            min_count: Minimum count for reaction edges (filters edges first)
            max_barrier: Maximum barrier height in kcal/mol

        Returns:
            Filtered pathways dictionary
        """
        valid_edges = set()
        if min_count is not None:
            for edge_key, edge in self.edge_map.items():
                start_id, end_id = map(int, edge_key.split("-"))
                if self.edge_meets_count_requirement(start_id, end_id, min_count):
                    valid_edges.add(edge_key)
        else:
            valid_edges = set(self.edge_map.keys())

        filtered = {}

        for node_id, pathways in self.pathways.items():
            valid_pathways = []

            # Sort by path length first, then by barrier height
            sorted_pathways = sorted(
                pathways,
                key=lambda p: (
                    p.get("path_length", 0),
                    p.get("barrier_height", float("inf")),
                ),
            )

            for pathway in sorted_pathways:
                if min_count is not None:
                    path = pathway["path"]
                    min_count_met = all(
                        self.edge_meets_count_requirement(
                            path[i], path[i + 1], min_count
                        )
                        for i in range(len(path) - 1)
                    )
                    if not min_count_met:
                        continue

                if (
                    max_barrier is not None
                    and pathway.get("barrier_height", float("inf")) > max_barrier
                ):
                    continue

                valid_pathways.append(pathway)

                if len(valid_pathways) >= 5:
                    break

            if valid_pathways:
                filtered[node_id] = valid_pathways

        return filtered

    def get_best_pathways(
        self, filtered_pathways: Optional[Dict] = None
    ) -> Dict[int, Dict]:
        """Get the best (lowest barrier) pathway to each node.

        Args:
            filtered_pathways: Pre-filtered pathways, or None to use all

        Returns:
            Dictionary mapping node_id to best pathway data
        """
        pathways_to_use = (
            filtered_pathways if filtered_pathways is not None else self.pathways
        )

        best_pathways = {}
        for node_id, pathways in pathways_to_use.items():
            if not pathways:
                continue

            best = min(pathways, key=lambda p: p["barrier_height"])
            best_pathways[node_id] = best

        return best_pathways

    def print_summary(self, filtered_pathways: Optional[Dict] = None):
        """Print a summary of pathway analysis."""
        pathways_to_use = (
            filtered_pathways if filtered_pathways is not None else self.pathways
        )

        total_nodes = len([p for p in pathways_to_use.values() if p])
        total_pathways = sum(len(p) for p in pathways_to_use.values())

        print("\nPathway Analysis Summary:")
        print(f"- Reachable nodes: {total_nodes}")
        print(f"- Total pathways: {total_pathways}")

        if total_nodes > 0:
            all_barriers = [
                p["barrier_height"]
                for pathways in pathways_to_use.values()
                for p in pathways
            ]
            if all_barriers:
                print(
                    f"- Barrier range: {min(all_barriers):.2f} - {max(all_barriers):.2f} kcal/mol"
                )
                print(f"- Average barrier: {np.mean(all_barriers):.2f} kcal/mol")

    def get_cached_profile(self, node_id: int, path_index: int = 0):
        """Get reaction profile with caching for repeated access."""
        cache_key = (node_id, path_index)

        if cache_key in self._profile_cache:
            return self._profile_cache[cache_key]

        if node_id not in self.pathways or path_index >= len(self.pathways[node_id]):
            return None

        profile = self.pathways[node_id][path_index]["profile"]
        self._profile_cache[cache_key] = profile

        return profile

    def get_edge(self, start_id: int, end_id: int) -> Optional[Dict]:
        """Get edge data from edge_map for fast lookup."""
        edge_key = f"{start_id}-{end_id}"
        return self.edge_map.get(edge_key)

    def get_edge_count(self, start_id: int, end_id: int) -> int:
        """Get edge count from edge_map for fast lookup."""
        edge = self.get_edge(start_id, end_id)
        return edge.get("count", 1) if edge else 1

    def get_edge_type(self, start_id: int, end_id: int) -> str:
        """Get edge type from edge_map for fast lookup."""
        edge = self.get_edge(start_id, end_id)
        return edge.get("type", "").lower() if edge else ""

    def edge_meets_count_requirement(
        self, start_id: int, end_id: int, min_count: int
    ) -> bool:
        """Check if an edge meets the count requirement.

        Only 'reaction' type edges need to meet the count threshold.
        All other edge types automatically pass regardless of count.

        Args:
            start_id: Starting node ID
            end_id: Ending node ID
            min_count: Minimum count requirement

        Returns:
            True if edge meets requirement, False otherwise
        """
        edge = self.get_edge(start_id, end_id)
        if not edge:
            return False

        edge_type = edge.get("type", "").lower()

        if edge_type != "reaction":
            return True

        # Reaction edges must meet the count threshold
        edge_count = edge.get("count", 1)
        return edge_count >= min_count

    def get_transition_state_data(self, start_id: int, end_id: int) -> Optional[Dict]:
        """Get transition state data for an edge."""
        edge = self.get_edge(start_id, end_id)
        if not edge or "ts" not in edge:
            return None

        try:
            ts_mol = Chem.MolFromMolBlock(edge["ts"])
            if not ts_mol:
                return None

            return {
                "mol": ts_mol,
                "mol_block": edge["ts"],
                "energy": edge.get("ts_energy", 0.0),
                "svg": self._generate_svg(self._prepare_drawing_mol(ts_mol)),
            }
        except Exception as e:
            print(f"Error processing TS data: {e}")
            return None

    def get_pathway_with_ts_details(
        self, node_id: int, path_index: int = 0
    ) -> Optional[Dict]:
        """Get detailed pathway information including all transition states."""
        if node_id not in self.pathways or path_index >= len(self.pathways[node_id]):
            return None

        pathway = self.pathways[node_id][path_index]
        path = pathway["path"]
        profile = pathway["profile"]

        detailed_steps = []

        for i, point in enumerate(profile):
            step_data = {
                "step": point["step"],
                "type": point["type"],
                "energy": point["energy"],
                "reaction_coordinate": point["reaction_coordinate"],
            }

            if point["type"] == "minimum":
                node_data = self.nodes[point["node_id"]]
                step_data.update(
                    {
                        "node_id": point["node_id"],
                        "svg": node_data.get("svg", ""),
                        "charge": node_data.get("charge", 0),
                        "weight": node_data.get("weight", 0),
                    }
                )

                if "edge_type" in point:
                    step_data["edge_type"] = point["edge_type"]

            elif point["type"] == "transition_state":
                step_data.update(
                    {"label": point["node_id"], "ts_data": point.get("ts_data", "")}
                )

                if "ts_data" in point:
                    try:
                        ts_mol = Chem.MolFromMolBlock(point["ts_data"])
                        if ts_mol:
                            step_data["svg"] = self._generate_svg(
                                self._prepare_drawing_mol(ts_mol)
                            )
                    except Exception:
                        pass

            detailed_steps.append(step_data)

        return {
            "path": path,
            "profile": profile,
            "detailed_steps": detailed_steps,
            "barrier_height": pathway["barrier_height"],
            "path_length": pathway["path_length"],
        }

    def to_visualization_data(self) -> Dict[str, Any]:
        """Convert ReactionNetwork data to format required by JavaScript
        visualization.

        This method provides the exact data structure expected by the existing
        JavaScript visualization code, ensuring compatibility while using the
        unified ReactionNetwork's advanced features.

        Returns:
            Dict containing data in JavaScript-compatible format
        """
        svg_map = {}
        for node_id, node_data in self.nodes.items():
            svg_map[node_id] = node_data.get("svg", "")

        mol_energy_map = {}
        for node_id, node_data in self.nodes.items():
            mol_energy_map[node_id] = node_data.get("energy", 0.0)

        charge_map = {}
        for node_id, node_data in self.nodes.items():
            charge_map[node_id] = node_data.get("charge", 0)

        weight_map = {}
        for node_id, node_data in self.nodes.items():
            weight_map[node_id] = node_data.get("weight", 0)

        rel_energy_map = {}
        substrate_energy = self.nodes[self.substrate_id]["energy"]
        for node_id, node_data in self.nodes.items():
            rel_energy = hartree2kcalmol(node_data["energy"] - substrate_energy)
            rel_energy_map[node_id] = rel_energy

        original_edge_data = []
        for edge in self.edges:
            edge_data = {
                "begin": edge["begin"],
                "end": edge["end"],
                "type": edge["type"],
                "count": edge.get("count", 0),
                "barrier": edge.get("barrier", 0.0),
                "rxn_energy": edge.get("rxn_energy", 0.0),
                "smaller_products": edge.get("smaller_products", []),
            }

            if "ts" in edge and edge["ts"]:
                edge_data["ts"] = edge["ts"]
                edge_data["ts_energy"] = edge.get("ts_energy")
                # TODO: here I need to also Thermoprocess the TS

            original_edge_data.append(edge_data)

        best_path_energy_map = {}
        best_path_barrier_map = {}
        best_path_map = {}

        for node_id, pathways in self.pathways.items():
            if pathways:
                best_pathway = pathways[0]
                best_path_energy_map[node_id] = best_pathway["profile"][-1][
                    "energy"
                ]  
                best_path_barrier_map[node_id] = best_pathway["barrier_height"]
                best_path_map[node_id] = best_pathway["path"]

        stereoisomer_groups = {}
        for node_id, node_data in self.nodes.items():
            if "sibling_ids" in node_data and node_data["sibling_ids"]:
                group_members = [node_id] + node_data["sibling_ids"]
                group_key = min(group_members)
                if group_key not in stereoisomer_groups:
                    stereoisomer_groups[group_key] = group_members

        pathway_lookup = {}
        pathway_profiles = {}
        pathway_barriers = {}

        for node_id, pathways in self.pathways.items():
            if pathways:
                pathway_lookup[node_id] = []
                pathway_profiles[node_id] = []
                pathway_barriers[node_id] = []

                for i, pathway in enumerate(pathways):
                    pathway_lookup[node_id].append(
                        {
                            "pathIndex": i,
                            "path": pathway["path"],
                            "pathLength": pathway["path_length"],
                            "barrierHeight": pathway["barrier_height"],
                        }
                    )

                    js_compatible_profile = []
                    for point in pathway["profile"]:
                        js_compatible_profile.append(
                            {
                                "nodeId": point["node_id"],
                                "energy": point["energy"],
                                "step": point["step"],
                                "reactionCoordinate": point.get(
                                    "reaction_coordinate", point["step"]
                                ),
                            }
                        )

                    pathway_profiles[node_id].append(js_compatible_profile)
                    pathway_barriers[node_id].append(pathway["barrier_height"])

        config_data = {
            "temperature": getattr(self.config, "temperature", 313.15),
            "R": 1.987e-3,  # Gas constant in kcal/(mol·K)
            "fitParams": getattr(self.config, "fit_params", {}),
            "pKaSlope": 0.43496246,
            "pKaIntercept": -114.66798804,
        }

        

        # Return data in format expected by JavaScript
        return {
            "svgMap": svg_map,
            "molEnergyMap": mol_energy_map,
            "chargeMap": charge_map,
            "weightMap": weight_map,
            "relEnergyMap": rel_energy_map,
            "originalEdgeData": original_edge_data,
            "substrateId": self.substrate_id,
            "energyType": self.energy_type,
            "temperature": config_data["temperature"],
            "pH": self.pH,
            "fitParams": config_data["fitParams"],
            "nodeLabels": self.node_labels,
            "bestPathEnergyMap": best_path_energy_map,
            "bestPathBarrierMap": best_path_barrier_map,
            "bestPathMap": best_path_map,
            "stereoisomerGroups": stereoisomer_groups,
            "precomputedPathways": self.pathways, 
            "pathwayLookup": pathway_lookup, 
            "pathwayProfiles": pathway_profiles, 
            "pathwayBarriers": pathway_barriers,  
            "configData": config_data,  
            "advancedFeatures": {
                "hasTransitionStates": any("ts" in edge for edge in self.edges),
                "hasReactionProfiles": len(self.pathways) > 0,
                "pathwayCount": sum(
                    len(pathways) for pathways in self.pathways.values()
                ),
                "reachableNodes": len([p for p in self.pathways.values() if p]),
                "hasSteReoisomers": len(stereoisomer_groups) > 0,
                "totalSteReoisomerGroups": len(stereoisomer_groups),
            },
        }

    def _load_nodes(self, level) -> None:
        """Load molecular nodes and compute their properties."""
        print("Loading and processing nodes...")

        node_files = list(self.config.node_data.glob("*.json"))
        nodes = [Node.load(f) for f in node_files]

        # Process template for alignment
        template_mol = None
        for n in nodes:
            if n.idx == self.substrate_id:
                mol = self._get_molecule_data(n, level)
                if mol:
                    template_mol = self._prepare_drawing_mol(mol)
                break

        for n in tqdm(nodes):
            mol = self._get_molecule_data(n, level)

            if not mol:
                print(f"Warning: Node {n.idx} has no molecular data")
                continue

            node_data = {
                "idx": n.idx,
                "mol": mol,
                "svg": self._generate_svg(self._prepare_drawing_mol(mol), template_mol),
                "energy": self._extract_energy(mol),
                "charge": Chem.rdmolops.GetFormalCharge(mol) if mol else 0,
                "weight": (
                    round(Chem.rdMolDescriptors.CalcExactMolWt(mol)) if mol else 0
                ),
            }
            xtb_mol = n.get_mol("gfn2-xtb")
            if xtb_mol.HasProp("sibling_labels"):
                slabels = eval(xtb_mol.GetProp("sibling_labels"))
                # convert labels to node indices
                sibling_ids = [
                    self.label2id[label] for label in slabels if label in self.label2id
                ]
                node_data["sibling_ids"] = sibling_ids
            node_data["origin-type"] = (
                xtb_mol.GetProp("origin-type") if xtb_mol.HasProp("origin-type") else ""
            )
            self.nodes[n.idx] = node_data

    def _get_molecule_data(self, node: Node, level: str) -> Optional[Chem.Mol]:
        """Get molecule data from a node."""
        standard_state = self.STANDARD_STATES.get(node.smi, 1.0)
        mol = node.get_qm(
            str(self.config.qm_data),
            level=level,
            T=float(self.config.temperature),
            c_M=standard_state,
        )
        if mol is None:
            mol = node.get_mol("gfn2-xtb")
        return mol

    @staticmethod
    def _prepare_drawing_mol(mol: Chem.Mol) -> Chem.Mol:
        """Prepare a molecule for drawing."""
        if mol is None:
            return None

        try:
            mol = Chem.RemoveHs(mol)
        except Exception as e:
            print(f"Error removing hydrogens: {e}")

        mol.RemoveAllConformers()
        return mol

    def _extract_energy(self, mol: Chem.Mol) -> float:
        """Extract energy from a molecule."""
        if mol is None:
            return 0.0

        try:
            return mol.GetDoubleProp(self.prop_name)
        except Exception as e:
            print(f"Error extracting energy from molecule: {e}")
            return 0.0

    def _generate_svg(
        self, mol: Chem.Mol, template_mol: Optional[Chem.Mol] = None
    ) -> str:
        """Generate SVG representation of a molecule."""
        if mol is None:
            return ""

        return Draw.MolsToGridImage(
            [mol], molsPerRow=1, subImgSize=(500, 400), useSVG=True
        )

    def _process_reactions(self) -> None:
        """Process reactions and create edge data."""
        print("Processing reactions...")

        reaction_files = list(self.config.reaction_data.glob("*.json"))
        reactions = [Reaction.load(f) for f in reaction_files]

        for reaction in tqdm(reactions):
            if reaction.rxn_type.lower() == "add":
                continue

            elif reaction.rxn_type.lower() == "reaction":
                if not reaction.reactant_ids:
                    print(f"Warning: Reaction {reaction.id} has no reactants")
                    continue

                reactant_id = reaction.reactant_ids[0]

                product_nodes = []
                for pid in reaction.product_ids:
                    node_files = list(self.config.node_data.glob(f"{pid}-*.json"))
                    if node_files:
                        product_nodes.append(Node.load(node_files[0]))

                if not product_nodes:
                    print(f"Warning: Reaction {reaction.id} has no valid product nodes")
                    continue

                product_nodes.sort(
                    key=lambda x: (
                        x.get_mol("gfn2-xtb").GetNumAtoms()
                        if x.get_mol("gfn2-xtb")
                        else 0
                    ),
                    reverse=True,
                )

                product_id = product_nodes[0].idx
                smaller_products = [n.idx for n in product_nodes[1:]]

                # Skip self-edges
                if reactant_id == product_id:
                    continue

                # Create edge
                edge = {
                    "id": len(self.edges),
                    "begin": reactant_id,
                    "end": product_id,
                    "type": reaction.rxn_type,
                    "count": (
                        reaction.count
                        if hasattr(reaction, "count") and reaction.count
                        else 1
                    ),
                }

                if smaller_products:
                    edge["smaller_products"] = smaller_products

                # See if there is a TS for this reaction
                ts_file = self.config.qm_data / f"ts-{reactant_id}-{product_id}.sdf"
                if ts_file.exists():
                    try:
                        with Chem.SDMolSupplier(
                            ts_file, removeHs=False, sanitize=False
                        ) as suppl:
                            if not suppl:
                                raise ValueError(
                                    f"No valid molecules found in {ts_file}"
                                )
                            # Assuming we want the first molecule in the SDF file
                            ts = next(suppl)

                        # Extract TS energy while we have the molecule with properties
                        if ts and ts.HasProp(self.prop_name):
                            edge["ts_energy"] = ts.GetDoubleProp(self.prop_name)

                        edge["ts"] = Chem.MolToMolBlock(ts)

                    except Exception as e:
                        print(
                            f"Error processing TS for {reactant_id}-{product_id}: {e}"
                        )

                self.edges.append(edge)
                edge_key = f"{reactant_id}-{product_id}"
                self.edge_map[edge_key] = edge

            else:
                # Handle other reaction types as direct edges
                if not reaction.reactant_ids or not reaction.product_ids:
                    print(
                        f"Warning: Reaction {reaction.id} missing reactants or products"
                    )
                    continue

                reactant_id = reaction.reactant_ids[0]
                product_id = reaction.product_ids[0]

                if reactant_id == product_id:
                    print(
                        f"⚠️  Skipping direct self-edge: {reactant_id} → {product_id}"
                    )
                    continue

                edge = {
                    "id": len(self.edges),
                    "begin": reactant_id,
                    "end": product_id,
                    "type": reaction.rxn_type,
                    "count": (
                        reaction.count
                        if hasattr(reaction, "count") and reaction.count
                        else 1
                    ),
                }
                self.edges.append(edge)
                edge_key = f"{reactant_id}-{product_id}"
                self.edge_map[edge_key] = edge

    def _calculate_relative_energies(self) -> None:
        """Calculate relative energies for all nodes based on absolute
        energies.

        All energies are relative to the substrate (substrate = 0).
        """
        print("Calculating relative node energies...")

        substrate_energy = self.nodes[self.substrate_id]["energy"]

        for node_id, node_data in self.nodes.items():
            rel_energy = hartree2kcalmol(node_data["energy"] - substrate_energy)
            node_data["rel_energy"] = rel_energy

        print(f"Calculated relative energies for {len(self.nodes)} nodes")

    def plot_reaction_profile(
        self, node_id: int, path_index: int = 0, show_ts: bool = True
    ):
        """Plot a reaction profile with optional transition states.

        Args:
            node_id: Target node ID
            path_index: Which path to plot (default: 0 = best path)
            show_ts: Whether to show transition states
        """
        try:
            import matplotlib.pyplot as plt
            from scipy.interpolate import CubicHermiteSpline
        except ImportError:
            print("matplotlib and scipy required for plotting")
            return None

        if node_id not in self.pathways or path_index >= len(self.pathways[node_id]):
            print(f"No pathway found for node {node_id}, path index {path_index}")
            return None

        pathway = self.pathways[node_id][path_index]
        profile = pathway["profile"]

        if not profile:
            print("Empty profile")
            return None

        coords = [p["reaction_coordinate"] for p in profile]
        energies = [p["energy"] for p in profile]

        # Separate minima and transition states
        minima_coords = [
            p["reaction_coordinate"] for p in profile if p["type"] == "minimum"
        ]
        minima_energies = [p["energy"] for p in profile if p["type"] == "minimum"]
        minima_labels = [str(p["node_id"]) for p in profile if p["type"] == "minimum"]

        ts_coords = [
            p["reaction_coordinate"] for p in profile if p["type"] == "transition_state"
        ]
        ts_energies = [p["energy"] for p in profile if p["type"] == "transition_state"]
        ts_labels = [p["node_id"] for p in profile if p["type"] == "transition_state"]

        # Create plot
        fig, ax = plt.subplots(figsize=(12, 6))

        # Interpolate smooth curve
        if len(coords) > 2:
            spline = CubicHermiteSpline(coords, energies, np.zeros_like(coords))
            fine_coords = np.linspace(min(coords) - 0.1, max(coords) + 0.1, 1000)
            fine_energies = spline(fine_coords)
            ax.plot(fine_coords, fine_energies, "b-", alpha=0.7, linewidth=2)

        # Plot minima
        ax.scatter(
            minima_coords, minima_energies, c="blue", s=80, zorder=5, label="Minima"
        )

        # Plot transition states
        if show_ts and ts_coords:
            ax.scatter(
                ts_coords,
                ts_energies,
                c="red",
                s=80,
                marker="x",
                zorder=5,
                label="Transition States",
            )

        # Add labels
        for i, (coord, energy, label) in enumerate(
            zip(minima_coords, minima_energies, minima_labels)
        ):
            ax.annotate(
                label,
                (coord, energy),
                xytext=(0, 10),
                textcoords="offset points",
                ha="center",
                fontsize=9,
            )

        if show_ts:
            for coord, energy, label in zip(ts_coords, ts_energies, ts_labels):
                ax.annotate(
                    label,
                    (coord, energy),
                    xytext=(0, 10),
                    textcoords="offset points",
                    ha="center",
                    fontsize=9,
                    color="red",
                )

        # Add barrier annotation
        barrier_height = pathway["barrier_height"]
        ax.text(
            0.02,
            0.98,
            f"Max Barrier: {barrier_height:.2f} kcal/mol",
            transform=ax.transAxes,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
        )

        ax.set_xlabel("Reaction Coordinate")
        ax.set_ylabel("Gibbs Free Energy (kcal/mol)")
        ax.set_title(f"Reaction Profile to Node {node_id}")
        ax.grid(True, alpha=0.3)
        ax.legend()

        plt.tight_layout()
        return fig, ax
