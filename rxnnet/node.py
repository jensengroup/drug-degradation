import copy
import hashlib
import json
from io import StringIO
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolHash
from rxnnet.thermo import Thermochemistry


def get_element_graph(mol):
    copy_mol = copy.deepcopy(mol)
    Chem.RemoveStereochemistry(copy_mol)
    m = Chem.RemoveAllHs(copy_mol)
    return rdMolHash.MolHash(m, rdMolHash.HashFunction.ElementGraph)


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
    def __init__(self, idx: int, data: dict = {}):
        self.idx = idx
        self._data = copy.deepcopy(data)  # {lot: Chem.Mol}

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
        }

    @classmethod
    def from_dict(cls, d):
        node = cls(d["idx"])
        for lot, sdf_str in d["data"].items():
            mol = mol_from_sdf_string(sdf_str)
            node._data[lot] = mol
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
