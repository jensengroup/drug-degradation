"""Utility functions for molecular processing and file management."""

import hashlib
import heapq
import json
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from tooltoad.chemutils import canonicalize_resonance, hartree2kcalmol

from d2.network import Node


def rm_db_stereo_info(mol):
    for bond in mol.GetBonds():
        if int(bond.GetStereo()) != 0:
            bond.SetStereo(Chem.rdchem.BondStereo.STEREOANY)


def standardize_mol(mol):
    Chem.SanitizeMol(mol)
    rm_db_stereo_info(mol)
    mol = canonicalize_resonance(mol)
    Chem.SanitizeMol(mol)
    return mol


def get_mol(atoms, coords, charge, remove_db_stereo=True):
    mol = Chem.MolFromXYZBlock(ac2xyz(atoms, coords))
    rdDetermineBonds.DetermineBonds(mol, charge=charge)
    mol = standardize_mol(mol)
    return mol


def calculate_mol_hash(mol: Chem.Mol) -> str:
    """Calculate molecular hash from SMILES string."""
    mol = standardize_mol(mol)
    smiles = Chem.MolToSmiles(mol)
    return hashlib.md5(smiles.encode()).hexdigest()[:12]


def make_mol(data: Dict[str, Any], coords: str = "opt_coords") -> Chem.Mol:
    """Create RDKit molecule from calculation data."""
    # Import here to avoid circular imports if needed
    from tooltoad.chemutils import ac2mol

    mol = ac2mol(data["atoms"], data[coords])
    rdDetermineBonds.DetermineBonds(mol, charge=data["charge"])

    mol = standardize_mol(mol)

    for prop in ["electronic_energy", "gibbs_energy"]:
        if prop in data:
            mol.SetDoubleProp(prop, data[prop])

    return mol


def ac2xyz(atoms: List[str], coords: np.ndarray) -> str:
    """Convert atoms and coordinates to XYZ format string."""
    xyz_lines = [str(len(atoms)), ""]  # Number of atoms and empty comment line
    for atom, coord in zip(atoms, coords):
        xyz_lines.append(f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")
    return "\n".join(xyz_lines)


def get_random_str() -> str:
    """Generate a random string of 6 characters."""
    return hashlib.md5(os.urandom(16)).hexdigest()[:6]


def track_nodes(track_file, node_files: List[Path]) -> List[Path]:
    """Track processed nodes to avoid reprocessing."""
    track_file = Path(track_file)
    processed_nodes = []

    if track_file.is_file():
        with open(track_file, "r") as f:
            track_data = f.read().strip().split("\n")
            processed_nodes = [
                line.split(",")[0] for line in track_data if line.strip()
            ]

        # Filter out already processed nodes
        node_files = [f for f in node_files if f.stem not in processed_nodes]
        print(f"Found {len(processed_nodes)} already processed nodes, skipping them.")
        print(f"Left with {len(node_files)} for processing")

    # Write processed nodes to the file
    with open(track_file, "a") as f:
        for file in node_files:
            f.write(f"{file.stem},\n")

    return node_files


def keep_strict(
    node_files: List[Path],
    reaction_data: Path,
    rxn_types: List[str] = ["protonation", "deprotonation", "tautomer"],
) -> List[Path]:
    """Filter nodes based on existing reactions for strict processing."""

    reaction_files = list(reaction_data.glob("*.json"))
    avoid = [
        f.stem.split("-")[1]
        for f in reaction_files
        if f.stem.split("-")[-1] in rxn_types
    ]

    node_files = [f for f in node_files if f.stem.split("-")[0] not in avoid]
    return node_files


class Reaction:
    """Represents a chemical reaction in the network."""

    def __init__(
        self,
        reactant_ids,
        product_ids,
        count=1,
        rxn_type: str = "reaction",
    ):
        self.reactant_ids = sorted(reactant_ids)
        self.product_ids = sorted(product_ids)
        self.id = f"{'_'.join([str(i) for i in self.reactant_ids])}-{'_'.join([str(i) for i in self.product_ids])}-{rxn_type.lower()}"
        self.count = count
        self.rxn_type = rxn_type

    def get_energy_profile(self):
        """Get the energy profile for this reaction."""

    def find_ts(self):
        """Find transition state for this reaction."""

    def run_irc(self):
        """Run IRC calculation for this reaction."""

    def plot(self):
        """Plot the reaction profile."""

    def show_ts(self):
        """Show the transition state structure."""

    def to_dict(self) -> Dict[str, Any]:
        """Convert reaction to dictionary."""
        return {
            "reactant_ids": self.reactant_ids,
            "product_ids": self.product_ids,
            "count": self.count,
            "rxn_type": self.rxn_type,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Reaction":
        """Create reaction from dictionary."""
        return cls(
            reactant_ids=data["reactant_ids"],
            product_ids=data["product_ids"],
            count=data["count"],
            rxn_type=data.get("rxn_type", "reaction"),
        )

    def save(self, directory: Path) -> str:
        """Save reaction to JSON file."""
        filename = f"{self.id}.json"
        with open(directory / filename, "w") as f:
            json.dump(self.to_dict(), f)
        return filename

    @classmethod
    def load(cls, filename: Path) -> "Reaction":
        """Load reaction from JSON file."""
        with open(filename, "r") as f:
            decoder = json.JSONDecoder()
            content = f.read()
            content = remove_multiline_json_keys(content, "json")
            data, _ = decoder.raw_decode(content)
        return cls.from_dict(data)


def remove_multiline_json_keys(text: str, key: str = "json") -> str:
    """Remove multiline JSON keys from text."""
    pattern = rf'"{re.escape(key)}"\s*:\s*(\{{|\[)'  # match "json": { or "json": [
    result = []
    i = 0

    while i < len(text):
        match = re.search(pattern, text[i:])
        if not match:
            result.append(text[i:])
            break

        start = i + match.start()
        open_char = match.group(1)
        close_char = "}" if open_char == "{" else "]"

        result.append(text[i:start])  # keep everything before "json":

        # Find the matching closing brace/bracket
        j = start + match.end() - match.start()
        depth = 1
        while j < len(text) and depth > 0:
            if text[j] == open_char:
                depth += 1
            elif text[j] == close_char:
                depth -= 1
            j += 1

        # Skip trailing whitespace and comma after the block
        while j < len(text) and text[j] in " \t\r\n":
            j += 1
        if j < len(text) and text[j] == ",":
            j += 1

        i = j  # continue after removed block

    return "".join(result)


def position_species(species_list: List[Any]) -> Tuple[List[str], np.ndarray]:
    """Position multiple species in 3D space to maximize separation."""
    mol_sizes = []
    all_atoms = []
    all_coords = []

    for species in species_list:
        mol = species.get_mol("gfn2-xtb")
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = np.asarray(mol.GetConformer().GetPositions())

        dm = np.zeros((len(atoms), len(atoms)))
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dm[i, j] = dm[j, i] = np.linalg.norm(coords[i] - coords[j])
        mol_sizes.append(np.max(dm))

        all_atoms.extend(atoms)
        all_coords.extend(coords)

    n_species = len(species_list)
    if n_species > 1:
        avg_size = np.mean(mol_sizes)
        radius = avg_size * (n_species ** (1 / 3))

        # Fibonacci sphere algorithm
        points = []
        phi = np.pi * (3.0 - np.sqrt(5.0))
        for i in range(n_species):
            y = 1 - (i / float(n_species - 1)) * 2
            radius_i = np.sqrt(1 - y * y)
            theta = phi * i
            x = np.cos(theta) * radius_i
            z = np.sin(theta) * radius_i
            points.append([x, y, z])

        for i in range(n_species):
            if i > 0:
                translation = np.array(points[i]) * radius
                start_idx = sum(len(s.atoms) for s in species_list[:i])
                end_idx = start_idx + len(species_list[i].atoms)
                for j in range(start_idx, end_idx):
                    all_coords[j] += translation

    return all_atoms, np.array(all_coords)


def fragment_charge(
    total_charge: int, n_fragments: int, min_max_frag_charge: Tuple[int, int] = (-2, 2)
):
    """Generate possible charge combinations for molecular fragments."""
    from itertools import product

    min_max_frag_charge = list(min_max_frag_charge)
    min_max_frag_charge[1] += 1
    min_max_frag_charge = range(*tuple(min_max_frag_charge))
    return [
        combo
        for combo in product(min_max_frag_charge, repeat=n_fragments)
        if sum(combo) == total_charge
    ]


def set_md_options(config):
    """Set up molecular dynamics options."""
    try:
        import click
        import typer
        from tooltoad.xtb import (
            MDOptions,
            MetaDynOptions,
            SCCOptions,
            WallOptions,
        )
    except ImportError:
        print("Required packages not available for MD options setup")
        return None

    options = [
        MDOptions(time=10, shake=0, temp=config.temperature),
        MetaDynOptions(kpush=0.15, alp=0.3),
        WallOptions(),
        SCCOptions(temp=9000),
    ]
    DEFAULT_OPTIONS = "\n".join([str(o) for o in options]) + "\n$cma"

    if typer.confirm(
        f"\nDo you want to change the MD options?\n\n{DEFAULT_OPTIONS}\n\n",
        default=False,
    ):
        PROMPT = "# Options for MD:"
        message = click.edit(PROMPT + "\n" + DEFAULT_OPTIONS)
        if message is not None:
            OPTIONS = "\n".join(
                [line for line in message.splitlines() if not line.startswith("#")]
            )
        else:
            OPTIONS = DEFAULT_OPTIONS
    else:
        OPTIONS = DEFAULT_OPTIONS

    return OPTIONS


def split_complex(
    atoms: List[str],
    coords: np.ndarray,
    overall_charge: int,
    xtb_options: Dict[str, Any] = None,
    n_cores: int = 1,
    scr: str = ".",
) -> Tuple[List[Chem.Mol], List[int], List[int]]:
    """Split a molecular complex into fragments and determine optimal
    charges."""
    import math

    from tooltoad.chemutils import ac2mol, hartree2kcalmol
    from tooltoad.xtb import xtb_calculate

    if xtb_options is None:
        xtb_options = {"alpb": "water"}

    complex_mol = ac2mol(atoms, coords, charge=overall_charge, use_xtb=True)
    frags = list(Chem.GetMolFrags(complex_mol, asMols=True))
    possible_charges = fragment_charge(overall_charge, len(frags))
    system_energies = []
    possible_multiplicities = []

    for pc in possible_charges:
        tmp_energies = []
        tmp_multiplicities = []
        for i, frag in enumerate(frags):
            frag_atoms = [a.GetSymbol() for a in frag.GetAtoms()]
            frag_coords = frag.GetConformer().GetPositions()
            charge = pc[i]
            multiplicity = (
                sum([a.GetAtomicNum() for a in frag.GetAtoms()]) - charge
            ) % 2 + 1
            tmp_multiplicities.append(multiplicity)

            results = xtb_calculate(
                frag_atoms,
                frag_coords,
                charge=charge,
                multiplicity=multiplicity,
                options=xtb_options,
                n_cores=n_cores,
                scr=scr,
            )

            if results["normal_termination"]:
                energy = results["electronic_energy"]
            else:
                energy = math.inf
            tmp_energies.append(energy)

        system_energies.append(sum(tmp_energies))
        possible_multiplicities.append(tmp_multiplicities)

    system_energies = np.asarray(system_energies)
    system_energies -= system_energies.min()
    system_energies *= hartree2kcalmol(1)
    best_charges = possible_charges[np.argmin(system_energies)]
    best_multiplicities = possible_multiplicities[np.argmin(system_energies)]

    return frags, best_charges, best_multiplicities


def choose_node(dir: Path, multiple: bool = False, types=["species", "complex"]):
    """Choose node files interactively."""
    from rich.columns import Columns
    from rich.console import Console
    from rich.panel import Panel
    from rich.prompt import Prompt

    console = Console()

    files = [f for f in dir.glob("*.json")]
    if not files:
        console.print("No nodes available in the network.")
        return
    species_files = sorted([f.name for f in files if int(f.name.split("-")[0]) < 1000])
    complex_files = sorted([f.name for f in files if int(f.name.split("-")[0]) > 999])

    # Create species and complex panels with files listed below them
    species_panel = Panel("\n".join(species_files), title="Species")
    complex_panel = Panel("\n".join(complex_files), title="Complexes")

    type2column = {"species": species_panel, "complex": complex_panel}

    # Create columns to display both panels side by side
    columns = Columns([type2column[t] for t in types])

    console.print(columns)
    if multiple:
        prompt = (
            "Select one or multiple species by index, separated by space (e.g. '1 2 3')"
        )
    else:
        prompt = "Select a species by index"
    user_input = Prompt.ask(prompt)
    if multiple:
        try:
            selected_ids = [int(idx) for idx in user_input.split(" ")]
        except ValueError:
            print(f"Invalid index provided: {user_input}")
            return
        filenames = []
        for idx in selected_ids:
            filename = [f for f in files if int(f.name.split("-")[0]) == idx]
            if len(filename) == 1:
                filenames.append(filename[0])
            else:
                print(f"Invalid index provided: {idx}")
        return filenames
    else:
        try:
            selected_id = int(user_input)
        except ValueError:
            print(f"Invalid index provided: {user_input}")
            return
        filename = [f for f in files if int(f.name.split("-")[0]) == selected_id]
        if len(filename) == 1:
            return filename[0]
        else:
            print(f"Invalid index provided: {user_input}")


def generate_and_save(
    generate_function,
    kwargs,
    origin: list[int],
    origin_type: str,
    new_nodes_dir,
    properties: dict = {},
    remove_db_stereo=True,
):
    """Generate new molecules and save them to the new_nodes directory."""
    results = generate_function(**kwargs)
    for result in results:
        if remove_db_stereo:
            rm_db_stereo_info(result)
            Chem.SanitizeMol(result)
        mol_hash = calculate_mol_hash(result)
        result.SetProp("origin-ids", str(sorted(origin)))
        result.SetProp("origin-type", origin_type)
        for name, value in properties.items():
            result.SetProp(name, str(value))
        with Chem.SDWriter(
            new_nodes_dir / f"{mol_hash}-{get_random_str()}.sdf"
        ) as writer:
            writer.write(result)
    return f"Generated {len(results)} molecules from {origin}."


# Energy filtering utilities


class NetworkEnergyCalculator:
    """Efficient calculator for relative energies in reaction networks."""

    def __init__(self, config):
        self.config = config
        self._mol_energy_map = None
        self._edges = None
        self._graph = None

    def _load_energy_data(self):
        """Load energy data once and cache it."""
        if self._mol_energy_map is not None:
            return

        self._mol_energy_map = {}

        # Load all nodes and their energies
        all_nodes = [Node.load(f) for f in self.config.node_data.glob("*.json")]
        for node in all_nodes:
            mol = node.get_qm(self.config.qm_data)
            if mol is None:
                mol = node.get_mol("gfn2-xtb")

            try:
                # Try to get the energy property (assuming gibbs energy)
                energy = mol.GetDoubleProp("l2l1_gibbs-energy") if mol else 0.0
            except (KeyError, AttributeError):
                try:
                    # Fallback to other energy properties
                    energy = mol.GetDoubleProp("energy") if mol else 0.0
                except (KeyError, AttributeError):
                    energy = 0.0

            self._mol_energy_map[node.idx] = energy

    def _build_network_graph(self):
        """Build the network graph once and cache it."""
        if self._graph is not None:
            return

        self._load_energy_data()

        reactions = [Reaction.load(f) for f in self.config.reaction_data.glob("*.json")]
        self._edges = []

        # Process reactions to create edges
        for reaction in reactions:
            if reaction.rxn_type.lower() == "add":
                continue
            elif reaction.rxn_type.lower() == "reaction":
                reactant_id = reaction.reactant_ids[0]
                product_nodes = [
                    Node.load(self.config.node_data.glob(f"{pid}-*.json").__next__())
                    for pid in reaction.product_ids
                ]
                product_nodes.sort(
                    key=lambda x: x.get_mol("gfn2-xtb").GetNumAtoms(), reverse=True
                )
                product_id = product_nodes[0].idx
                smaller_products = [n.idx for n in product_nodes[1:]]

                # Calculate deltaE
                delta_e = hartree2kcalmol(
                    self._mol_energy_map.get(product_id, 0)
                    + sum([self._mol_energy_map.get(pf, 0) for pf in smaller_products])
                    - self._mol_energy_map.get(reactant_id, 0)
                    + self.config.offsets.get(
                        self.config.solvent.lower() if self.config.solvent else "", {}
                    ).get(reaction.rxn_type.lower(), 0)
                )

                edge = {
                    "begin": reactant_id,
                    "end": product_id,
                    "deltaE": delta_e,
                }
                self._edges.append(edge)
            else:
                # Direct edge for single product reactions
                reactant_id = reaction.reactant_ids[0]
                product_id = reaction.product_ids[0]

                delta_e = hartree2kcalmol(
                    self._mol_energy_map.get(product_id, 0)
                    - self._mol_energy_map.get(reactant_id, 0)
                    + self.config.offsets.get(
                        self.config.solvent.lower() if self.config.solvent else "", {}
                    ).get(reaction.rxn_type.lower(), 0)
                )

                edge = {
                    "begin": reactant_id,
                    "end": product_id,
                    "deltaE": delta_e,
                }
                self._edges.append(edge)

        # Build adjacency list with edge weights (deltaE values)
        self._graph = defaultdict(list)
        for edge in self._edges:
            self._graph[edge["begin"]].append((edge["end"], edge["deltaE"]))
            self._graph[edge["end"]].append((edge["begin"], -edge["deltaE"]))

    def calculate_relative_energies(self, substrate_id: int) -> Dict[int, float]:
        """Calculate relative energies from substrate using shortest path
        algorithm."""
        self._build_network_graph()

        # Use Dijkstra's algorithm to find shortest energy paths from substrate
        distances = {substrate_id: 0.0}
        visited = set()
        queue = [(0.0, substrate_id)]  # (distance, node_id)

        while queue:
            current_dist, current_node = heapq.heappop(queue)

            if current_node in visited:
                continue

            visited.add(current_node)

            for neighbor, edge_weight in self._graph[current_node]:
                if neighbor not in visited:
                    new_dist = current_dist + edge_weight

                    if neighbor not in distances or new_dist < distances[neighbor]:
                        distances[neighbor] = new_dist
                        heapq.heappush(queue, (new_dist, neighbor))

        return distances


def filter_nodes_by_energy(
    node_files: List[Path], config, substrate_id: int, energy_threshold: float
) -> List[Path]:
    """Filter node files by relative energy threshold."""
    if energy_threshold == float("inf"):
        return node_files

    calculator = NetworkEnergyCalculator(config)
    relative_energies = calculator.calculate_relative_energies(substrate_id)

    # Filter nodes based on energy threshold
    filtered_files = []
    print(
        f"Filtering {len(node_files)} nodes with energy threshold: {energy_threshold} kcal/mol"
    )
    for file in node_files:
        node = Node.load(file)
        rel_energy = relative_energies.get(node.idx, float("inf"))
        if rel_energy <= energy_threshold:
            filtered_files.append(file)

    print(f"Filtered to {len(filtered_files)} nodes within energy threshold")
    return filtered_files


def select_and_filter_nodes(
    config,
    node_ids: Optional[Union[List[int], str]] = None,
    energy_threshold: Optional[float] = None,
    substrate_id: int = 1,
    types: Optional[List[str]] = None,
    track_file: Optional[str] = None,
    additional_filtering: Optional[Callable] = None,
) -> List[Path]:
    """Unified node selection and filtering logic for CLI commands."""
    if types is None:
        types = ["species"]

    # Step 1: Initial node selection
    if node_ids is None:
        # Interactive selection
        node_files = [
            Path(f) for f in choose_node(config.node_data, multiple=True, types=types)
        ]
        if not node_files:
            return []
    elif isinstance(node_ids, list) and not node_ids:
        node_files = list(config.node_data.glob("*.json"))
        print(f"Found {len(node_files)} total nodes.")
    elif isinstance(node_ids, str) and (
        node_ids.lower() == "all" or node_ids.strip() == ""
    ):
        node_files = list(config.node_data.glob("*.json"))
        print(f"Found {len(node_files)} total nodes.")
    elif isinstance(node_ids, list):
        node_files = []
        for node_id in node_ids:
            matching_files = list(config.node_data.glob(f"{node_id}-*.json"))
            if matching_files:
                node_files.append(matching_files[0])
            else:
                print(f"Warning: No file found for node ID {node_id}")
    else:
        try:
            ids = [
                int(x.strip()) for x in node_ids.replace(",", " ").split() if x.strip()
            ]
            node_files = []
            for node_id in ids:
                matching_files = list(config.node_data.glob(f"{node_id}-*.json"))
                if matching_files:
                    node_files.append(matching_files[0])
                else:
                    print(f"Warning: No file found for node ID {node_id}")
        except ValueError:
            print(f"Invalid node ID format: {node_ids}")
            return []

    if not node_files:
        print("No valid node files found.")
        return []

    # Apply additional command-specific filtering
    if additional_filtering:
        node_files = additional_filtering(node_files)

    # Apply tracking to avoid reprocessing
    if track_file:
        original_count = len(node_files)
        node_files = track_nodes(track_file, node_files)
        if len(node_files) < original_count:
            print(
                f"Skipped {original_count - len(node_files)} already processed nodes."
            )

    if not node_files:
        print("No nodes remaining after filtering and tracking.")
        return []

    print(f"Final selection: {len(node_files)} nodes to process.")
    return node_files
