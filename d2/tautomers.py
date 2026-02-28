import copy
import os

from rdkit import Chem
from rdkit.Chem import (
    ResonanceMolSupplier,
    rdDetermineBonds,
    rdDistGeom,
    rdMolHash,
    rdmolops,
)
from tooltoad.chemutils import (
    ac2mol,
    ac2xyz,
    canonicalize_solvent,
    read_multi_xyz,
)
from tooltoad.utils import WorkingDir, check_executable, stream
from tooltoad.xtb import xtb_calculate


def get_element_graph(mol):
    copy_mol = copy.deepcopy(mol)
    Chem.RemoveStereochemistry(copy_mol)
    m = Chem.RemoveAllHs(copy_mol)
    return rdMolHash.MolHash(m, rdMolHash.HashFunction.ElementGraph)


class TautomerGenerator:
    def __init__(self, crest_executable: str = "crest"):
        if check_executable(crest_executable) is not None:
            raise RuntimeError("CREST executable not found. Please install CREST.")
        self.crest_executable = crest_executable

    def __call__(
        self,
        mol: Chem.Mol,
        solvent: str | None = None,
        n_cores: int = 1,
    ):
        scr = os.getenv("SCRATCH", ".")
        mol = Chem.AddHs(mol)
        charge = rdmolops.GetFormalCharge(mol)
        rdDistGeom.EmbedMolecule(mol, randomSeed=42)
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = mol.GetConformer().GetPositions()
        opt_options = {"opt": None}
        if solvent is not None:
            opt_options["alpb"] = canonicalize_solvent(solvent, qm="xtb")
        opt_results = xtb_calculate(atoms, coords, charge=charge, options=opt_options)
        if not opt_results["normal_termination"]:
            raise RuntimeError(f"XTB optimization failed:\n{opt_results['log']}")
        opt_coords = opt_results["opt_coords"]
        work_dir = WorkingDir(root=scr)
        with open(work_dir / "input.xyz", "w") as f:
            f.write(ac2xyz(atoms, opt_coords))

        cmd = f"{self.crest_executable} input.xyz --chrg {charge} --T {n_cores} --tautomerize"
        if solvent is not None:
            cmd += f" --alpb {canonicalize_solvent(solvent, qm='xtb')}"
        cmd += " | tee tautomers.log"
        lines = list(stream(cmd, cwd=str(work_dir)))
        normal_termination = False
        for line in reversed(lines):
            if "CREST terminated normally" in line:
                normal_termination = True
                break
        if not normal_termination:
            raise RuntimeError(f"CREST failed to run:\n{lines}")
        try:
            atoms, coords, energies = read_multi_xyz(
                work_dir / "tautomers.xyz",
                extract_property_function=lambda x: float(x),
            )
            input_mol = ac2mol(atoms[0], opt_coords)
            query_mol = Chem.MolFromSmarts(Chem.MolToSmarts(Chem.RemoveHs(input_mol)))
        except Exception as e:
            print(e)
        finally:
            work_dir.cleanup()
        out_mols = [ac2mol(a, c) for a, c in zip(atoms, coords)]
        for om, energy in zip(out_mols, energies):
            rdDetermineBonds.DetermineBonds(om, charge=charge)
            om.SetDoubleProp("electronic_energy", energy)
        tautomers = [
            om for om in out_mols if get_element_graph(om) == get_element_graph(mol)
        ]
        # sort based on electronic energy
        tautomers.sort(key=lambda m: m.GetDoubleProp("electronic_energy"))
        # only keep tautomers with unique smiles
        unique_smiles = set()
        tautomers = [
            t
            for t in tautomers
            if Chem.MolToSmiles(t) not in unique_smiles
            and not unique_smiles.add(Chem.MolToSmiles(t))
        ]
        # get canonical resonance forms
        tautomers = [ResonanceMolSupplier(mol).__next__() for mol in tautomers]

        clean_mols = []
        for mol in tautomers:
            tmp_mol = copy.deepcopy(mol)
            for b in tmp_mol.GetBonds():
                b.SetBondType(Chem.BondType.SINGLE)
            degree_changes = []
            for input_idx, mol_idx in zip(
                input_mol.GetSubstructMatch(query_mol),
                tmp_mol.GetSubstructMatch(query_mol),
            ):
                if (
                    input_mol.GetAtomWithIdx(input_idx).GetDegree()
                    != tmp_mol.GetAtomWithIdx(mol_idx).GetDegree()
                ):
                    degree_changes.append(
                        (
                            input_mol.GetAtomWithIdx(input_idx).GetDegree(),
                            tmp_mol.GetAtomWithIdx(mol_idx).GetDegree(),
                        )
                    )
            if len(degree_changes) <= 2:
                clean_mols.append(mol)

        return clean_mols
