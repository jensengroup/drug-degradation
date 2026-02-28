import copy
import hashlib
import json
import re
from io import StringIO
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolHash, rdmolops
from tooltoad.chemutils import ac2mol
from tooltoad.thermo import Thermochemistry


def get_element_graph(mol):
    copy_mol = copy.deepcopy(mol)
    Chem.RemoveStereochemistry(copy_mol)
    m = Chem.RemoveAllHs(copy_mol)
    return rdMolHash.MolHash(m, rdMolHash.HashFunction.ElementGraph)


class Ensemble:
    def __init__(self, mols: list[Chem.Mol]):
        self.mols = mols

    @staticmethod
    def same_element_graph(mols: list[Chem.Mol], smi: str | None = None):
        if smi is None:
            test_smi = get_element_graph(mols[0])
        else:
            test_smi = get_element_graph(Chem.MolFromSmiles(smi))
        mask = []
        for mol in mols:
            mask.append(get_element_graph(mol) == test_smi)
        return mask

    @classmethod
    def from_xyz(
        cls,
        xyz_file,
        charge: int = 0,
        property_function: exec = lambda x: float(x),
        level: str = "energy",
        smi: str | None = None,
        n_skip: int = 0,
    ):
        with open(xyz_file, "r") as f:
            lines = f.readlines()

        n_atoms = int(lines[0])
        frame_size = n_atoms + 2 + n_skip  # Atoms + comment + atom count line
        n_frames = len(lines) // frame_size

        atoms = []
        coords = []
        property = []
        for i in range(n_frames):
            start = i * frame_size + 2  # Skip atom count and metadata lines

            property.append(property_function(lines[start - 1]))
            frame_coords = [
                list(map(float, line.split()[1:]))
                for line in lines[start : start + n_atoms]
            ]
            coords.append(frame_coords)
            atoms.append([line.split()[0] for line in lines[start : start + n_atoms]])
        for a, c, e in zip(atoms, coords, property):
            # check that it would get the same smiles
            mol = ac2mol(a, c)
            rdDetermineBonds.DetermineBonds(mol, charge=charge)

            conf = mol.GetConformer()

            conf.SetDoubleProp(level, e)
            mol.AddConformer(conf, assignId=True)

        return cls(mol)

    def to_sdf(self, file_path: str):
        with Chem.SDWriter(file_path) as writer:
            for mol in self.mols:
                writer.write(mol)

    @classmethod
    def from_sdf(cls, file_path: str):
        suppl = Chem.SDMolSupplier(file_path, removeHs=False)
        mols = [mol for mol in suppl if mol is not None]
        if not mols:
            raise ValueError(f"No valid molecules found in {file_path}")
        return cls(mols)

    def to_ac(self):
        atoms, coords = [], []
        for mol in self.mols:
            atoms.append([a.GetSymbol() for a in mol.GetAtoms()])
            coords.append(mol.GetConformer().GetPositions().tolist())
        return atoms, coords


def ac2id(atoms, coords, stero=False):
    tmp = ac2mol(atoms, coords, use_xtb=True)
    if stero:
        rdmolops.AssignStereochemistryFrom3D(tmp)
    smi = Chem.MolToSmiles(tmp)
    return hashlib.md5(smi.encode()).hexdigest()[:12]


LOT = "gfn2-xtb"


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
            assert (
                "Single_Point_Data" in self.l2["Geometries"][-1]
            ), "L2 file missing Single_Point_Data"
            assert self.l2["Geometries"][-1]["Single_Point_Data"][
                "Converged"
            ], "L2 file Single_Point_Data not converged"
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


def mol_to_sdf_string(mol: Chem.Mol) -> str:
    buf = StringIO()
    writer = Chem.SDWriter(buf)
    writer.write(mol)
    writer.flush()
    return buf.getvalue()


def mol_from_sdf_string(sdf_string: str) -> Chem.Mol:
    ms = Chem.SDMolSupplier()
    ms.SetData(sdf_string, removeHs=False)
    mol = next(ms)
    return mol


class Node:
    def __init__(self, idx: int, data: dict = {}, conformer_ensemble=None):
        self.idx = idx
        self._data = copy.deepcopy(data)  # {lot: Chem.Mol}
        self.conformer_ensemble = conformer_ensemble
        self.processes = []

    def __repr__(self):
        return f"<Node idx={self.idx}, lots={self.list_lots()}>"

    @property
    def smi(self):
        if not self._data:
            raise ValueError("No molecules available to get SMILES.")
        # Return the SMILES of the first molecule
        return Chem.MolToSmiles(next(iter(self._data.values())))

    @property
    def id(self):
        return hashlib.md5(self.smi.encode()).hexdigest()[:12]

    def add_mol(self, lot: str, mol: Chem.Mol):
        if not isinstance(mol, Chem.Mol):
            raise TypeError("Data must be a Chem.Mol.")
        self._data[lot] = mol

    def get_mol(self, lot: str) -> Chem.Mol:
        opt_lot = lot.split("/")[-1]
        for lvl in [lot, opt_lot]:
            try:
                return self._data[lvl]
            except KeyError:
                continue
        raise KeyError(f"No data found for level of theory '{lot}'")

    def get_qm(
        self, qm_data_path: str, level: str = "quick", T: float = 298.15, c_M: float = 1
    ) -> Chem.Mol:
        qm_data_path = Path(qm_data_path)
        f = qm_data_path / f"{self.idx}-{self.id}-{level}.sdf"
        if not f.is_file():
            return None
        with Chem.SDMolSupplier(qm_data_path / f, removeHs=False) as suppl:
            if not suppl:
                raise ValueError(f"No valid molecules found in {qm_data_path}")
            mol = next(suppl)
        # now get the qm data
        f_l1 = qm_data_path / f"{self.idx}-{self.id}-freq-{level}.json"
        f_l2 = qm_data_path / f"{self.idx}-{self.id}-sp-{level}.json"
        if not f_l1.is_file() or not f_l2.is_file():
            return None
        print(
            f"Calculating Gibbs energy for node {self.idx} at level {level} with T={T} and c_M={c_M}"
        )
        result = QMResult.from_files(f_l1, f_l2)
        gibbs_energy = result.get_energy(T=T, c_M=c_M, p_atm=None)
        for prop_name in mol.GetPropsAsDict().keys():
            if "energy" in prop_name.lower():
                print("Removing old energy property:", prop_name)
                mol.ClearProp(prop_name)
        mol.SetDoubleProp("l2l1_gibbs-energy", gibbs_energy)
        return mol

    def has_lot(self, lot: str) -> bool:
        return lot in self._data

    def list_lots(self) -> list[str]:
        return list(self._data.keys())

    def to_dict(self):
        return {
            "idx": self.idx,
            "data": {lot: mol_to_sdf_string(mol) for lot, mol in self._data.items()},
            "conformer_ensemble": (
                self.conformer_ensemble.to_dict() if self.conformer_ensemble else None
            ),
        }

    @classmethod
    def from_dict(cls, d):
        node = cls(d["idx"])
        for lot, sdf_str in d["data"].items():
            mol = mol_from_sdf_string(sdf_str)
            node.add_mol(lot, mol)
        if d.get("conformer_ensemble"):
            node.conformer_ensemble = Ensemble.from_dict(d["conformer_ensemble"])
        return node

    def save(self, dir: str | Path):
        filename = Path(f"{self.idx}-{self.id}.json")
        with open(dir / filename, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: Path):
        with open(path, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)


def remove_multiline_json_keys(text, key="json"):
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
