import json
import logging
import os
import tempfile
from pathlib import Path
from typing import List
from typing import Any, Dict, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem  # noqa: F401
from rdkit.Chem import (
    ResonanceMolSupplier,
    rdDetermineBonds,
    rdmolops,
)
from rdkit.Chem.rdchem import Mol
import hashlib
from itertools import product


from rxnnet.utils import stream, get_random_str

logger = logging.getLogger(__name__)


VDW_RADII = {"C": 1.7, "N": 1.55, "O": 1.52, "H": 1.2, "S": 1.8, "F": 1.47, "Cl": 1.75}

# convert solvent names to canonical names
# at least in xtb 6.7.1 the DCM name doesn't work
CANONICAL_SOLVENT_NAMES = {"xtb": {"dcm": "ch2cl2"}, "orca": {}}


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


def get_mol(atoms, coords, charge, remove_db_stereo=True, strict=True):
    """Build an RDKit Mol from atoms/coords with bond perception.

    If *strict* is True (default), raises on failure. If *strict* is
    False, returns None when bond perception or sanitization fails.
    """
    try:
        mol = Chem.MolFromXYZBlock(ac2xyz(atoms, coords))
        if mol is None:
            if strict:
                raise ValueError("MolFromXYZBlock returned None")
            return None
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
        Chem.SanitizeMol(mol)
        if remove_db_stereo:
            rm_db_stereo_info(mol)
            Chem.SanitizeMol(mol)
        return mol
    except Exception:
        if strict:
            raise
        return None


def get_multiplicity(mol, charge=None):
    """Return the lowest valid spin multiplicity for a molecule and charge.

    Singlet (1) when the electron count is even, doublet (2) when odd.
    *charge* overrides the formal charge read from *mol*.
    """
    if charge is None:
        charge = Chem.rdmolops.GetFormalCharge(mol)
    n_electrons = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge
    return 1 if n_electrons % 2 == 0 else 2


def calculate_mol_hash(mol: Chem.Mol) -> str:
    """Calculate molecular hash from SMILES string."""
    smiles = get_smiles(mol)
    return hashlib.md5(smiles.encode()).hexdigest()[:12]


def get_smiles(mol):
    mol = standardize_mol(mol)
    return Chem.MolToSmiles(mol)


def make_mol(data: Dict[str, Any], coords: str = "opt_coords") -> Chem.Mol:
    """Create RDKit molecule from calculation data."""
    # Import here to avoid circular imports if needed

    mol = ac2mol(data["atoms"], data[coords], perceive_connectivity=False)
    rdDetermineBonds.DetermineBonds(mol, charge=data["charge"])

    mol = standardize_mol(mol)

    for prop in ["electronic_energy", "gibbs_energy"]:
        if prop in data:
            mol.SetDoubleProp(prop, data[prop])

    return mol


def fragment_charge(
    total_charge: int, n_fragments: int, min_max_frag_charge: Tuple[int, int] = (-2, 2)
):
    """Generate possible charge combinations for molecular fragments."""

    min_max_frag_charge = list(min_max_frag_charge)
    min_max_frag_charge[1] += 1
    min_max_frag_charge = range(*tuple(min_max_frag_charge))
    return [
        combo
        for combo in product(min_max_frag_charge, repeat=n_fragments)
        if sum(combo) == total_charge
    ]


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


def canonicalize_resonance(mol):
    try:
        mol = ResonanceMolSupplier(mol).__next__()
    except StopIteration:
        pass
    return mol


def canonicalize_solvent(solvent: str, qm: str):
    assert qm.lower() in ["xtb", "orca"], "QM must be either xtb or orca"
    if solvent:
        return CANONICAL_SOLVENT_NAMES[qm.lower()].get(solvent.lower(), solvent)


class Constraint:
    def __init__(self, ids: list[int], value: float):
        self.ids = ids
        self.value = value
        assert all(isinstance(i, int) for i in ids)

    def __repr__(self):
        return f"{self.xtb_type.capitalize()} Constraint at {self.xtb_value} for ({','.join([str(x) for x in self.ids])})"

    @property
    def xtb_type(self):
        if len(self.ids) == 2:
            return "distance"
        elif len(self.ids) == 3:
            return "angle"
        elif len(self.ids) == 4:
            return "dihedral"
        else:
            raise ValueError

    @property
    def xtb_ids(self):
        return [i + 1 for i in self.ids]

    @property
    def xtb_value(self):
        return self.value if self.value else "auto"

    @property
    def xtb(self):
        return f"{self.xtb_type}: {', '.join([str(x) for x in self.xtb_ids])}, {self.xtb_value}"

    @property
    def orca_type(self):
        if len(self.ids) == 2:
            return "B"
        elif len(self.ids) == 3:
            return "A"
        elif len(self.ids) == 4:
            return "D"
        else:
            raise ValueError

    @property
    def orca_ids(self):
        return self.ids

    @property
    def orca_value(self):
        return self.value if self.value else ""

    @property
    def orca(self):
        return f"{{ {self.orca_type} {' '.join([str(x) for x in self.orca_ids])} {self.orca_value} C }}"


def same_connectivity(
    mol: Mol,
    atoms: list[str],
    opt_coords: list[list[float]],
    charge: int,
    multiplicity: int,
    scr: str,
) -> tuple[bool, None | np.ndarray]:
    """Check if the connectivity of the molecule is the same before and after
    optimization."""
    ac1 = rdmolops.GetAdjacencyMatrix(mol)
    new_mol = ac2mol(atoms, opt_coords, charge, multiplicity, scr, sanitize=False)
    ac2 = rdmolops.GetAdjacencyMatrix(new_mol)
    if (ac1 == ac2).all():
        return True, None
    else:
        logger.debug("Connectivity changed")
        return False, ac2 - ac1


def hartree2kcalmol(hartree: float) -> float:
    """Converts Hartree to kcal/mol."""
    return hartree * 627.509474


def read_multi_xyz(
    xyz_traj_file: str, extract_property_function: None | (str) = None, n_skip: int = 0
) -> tuple:
    """Reads a multi-frame XYZ trajectory file and returns a list of
    coordinates and optionally properties."""
    with open(xyz_traj_file, "r") as f:
        lines = f.readlines()

    n_atoms = int(lines[0])
    frame_size = n_atoms + 2 + n_skip  # Atoms + comment + atom count line
    n_frames = len(lines) // frame_size

    atoms = []
    coords = []
    property = []
    for i in range(n_frames):
        start = i * frame_size + 2  # Skip atom count and metadata lines
        if extract_property_function:
            property.append(extract_property_function(lines[start - 1]))
        frame_coords = [
            list(map(float, line.split()[1:]))
            for line in lines[start : start + n_atoms]
        ]
        coords.append(frame_coords)
        atoms.append([line.split()[0] for line in lines[start : start + n_atoms]])
    if extract_property_function:
        return atoms, coords, property
    return atoms, coords


# @hide_warnings
def _determineConnectivity(mol, usextb=False, **kwargs):
    """Determine bonds in molecule."""
    if usextb:
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = mol.GetConformer().GetPositions()
        charge = rdmolops.GetFormalCharge(mol)
        multiplicity = kwargs.get("multiplicity", 1)
        scr = kwargs.get("scr", ".")
        adj = gfnff_connectivity(atoms, coords, charge, multiplicity, scr)
        emol = Chem.EditableMol(mol)
        for i, j in np.argwhere(adj):
            if i > j:
                emol.AddBond(int(i), int(j), Chem.BondType.SINGLE)
        mol = emol.GetMol()
    else:
        try:
            rdDetermineBonds.DetermineConnectivity(mol, **kwargs)
        finally:
            # cleanup extended hueckel files
            try:
                os.remove("nul")
                os.remove("run.out")
            except FileNotFoundError:
                pass
    return mol


def xyz2ac(xyzblock: str):
    """Converts atom symbols and coordinates to xyz string."""
    lines = xyzblock.split("\n")
    atoms = []
    coords = []
    for line in lines[2:]:
        line = line.strip()
        if len(line) > 0:
            atom, x, y, z = line.split()
            atoms.append(atom)
            coords.append([float(x), float(y), float(z)])
        else:
            break
    return atoms, coords


def ac2xyz(atoms: List[str], coords: List[list]):
    """Converts atom symbols and coordinates to xyz string."""
    xyz = f"{len(atoms)}\n\n"
    for atom, coord in zip(atoms, coords):
        xyz += f"{atom} {coord[0]:.8f} {coord[1]:.8f} {coord[2]:.8f}\n"
    return xyz


def ac2mol(
    atoms: List[str],
    coords: List[list],
    charge: int = 0,
    multiplicity: int = 1,
    scr: str = ".",
    perceive_connectivity: bool = True,
    use_xtb: bool = True,
    sanitize: bool = False,
):
    """Converts atom symbols and coordinates to RDKit molecule."""
    xyz = ac2xyz(atoms, coords)
    rdkit_mol = Chem.MolFromXYZBlock(xyz)
    if sanitize:
        Chem.SanitizeMol(rdkit_mol)
    if perceive_connectivity:
        rdkit_mol = _determineConnectivity(
            rdkit_mol, usextb=use_xtb, charge=charge, multiplicity=multiplicity, scr=scr
        )
    return rdkit_mol


def gfnff_connectivity(atoms, coords, charge, multiplicity, scr):
    # Determine connectivity based on GFNFF-xTB implementation
    calc_dir = tempfile.TemporaryDirectory(dir=scr)
    tmp_file = Path(calc_dir.name) / "input.xyz"
    with open(tmp_file, "w") as f:
        f.write(ac2xyz(atoms, coords))
    CMD = f"xtb --gfnff {str(tmp_file.name)} --chrg {charge} --uhf {multiplicity - 1} --norestart --wrtopo blist"
    _ = list(stream(CMD, cwd=calc_dir.name))
    with open(Path(calc_dir.name) / "gfnff_lists.json", "r") as f:
        data_dict = json.load(f)
    calc_dir.cleanup()
    blist = data_dict["blist"]
    adj = np.zeros((len(atoms), len(atoms)), dtype=int)
    for i, j, _ in blist:
        adj[i - 1, j - 1] = 1
        adj[j - 1, i - 1] = 1
    return adj
