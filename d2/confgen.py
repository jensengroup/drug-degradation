import copy

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdDistGeom, rdForceFieldHelpers
from tooltoad.chemutils import (
    ConformerCalculator,
    ac2xyz,
    canonicalize_solvent,
    energy_filter_conformer,
)
from tooltoad.orca import orca_calculate
from tooltoad.utils import check_executable
from tooltoad.xtb import xtb_calculate


class ConformerGenerator:
    def __init__(self, mode: str, crest_executable: str = "crest"):
        self.mode = mode.lower()
        if self.mode not in ["etkdg", "crest"]:
            raise ValueError(f"Invalid mode: {self.mode}. Must be 'etkdg' or 'crest'.")
        if self.mode == "crest" and check_executable(crest_executable) is not None:
            raise RuntimeError("CREST executable not found. Please install CREST.")
        self.crest_executable = crest_executable

    def __call__(
        self,
        mol: Chem.Mol,
        solvent: str | None = None,
        n_cores: int = 1,
        scr: str = ".",
        randomSeed: int = 42,
        numConfs: int = 100,
        enforce_stereo: bool = True,
    ):
        assert (
            self.mode == "etkdg"
        ), "Only ETKDG mode is currently implemented in this class."
        mol = Chem.AddHs(mol)
        original_stereo = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if "?" in [stereo for _, stereo in original_stereo]:
            enforce_stereo = False
        original_smi = Chem.MolToSmiles(mol)
        rdDistGeom.EmbedMultipleConfs(
            mol,
            randomSeed=randomSeed,
            numConfs=numConfs,
            numThreads=n_cores,
            pruneRmsThresh=0.2,
        )
        if rdForceFieldHelpers.MMFFHasAllMoleculeParams(mol):
            rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(
                mol, numThreads=n_cores, maxIters=1000
            )
        elif rdForceFieldHelpers.UFFHasAllMoleculeParams(mol):
            rdForceFieldHelpers.UFFOptimizeMoleculeConfs(
                mol, numThreads=n_cores, maxIters=1000
            )
        else:
            raise RuntimeError(
                "No force field parameters available for the molecule. "
                "Please ensure the molecule has valid MMFF or UFF parameters."
            )
        opt_options = {"opt": None}
        if solvent is not None:
            xtb_solvent = canonicalize_solvent(solvent, "xtb")
            opt_options["alpb"] = xtb_solvent

        xtb_opt = ConformerCalculator(xtb_calculate, options=opt_options, scr=scr)
        mol, _ = xtb_opt(mol, n_cores=n_cores, memory=n_cores * 2)
        #  check that only the correct stereochemistry is present
        # filter
        mol = energy_filter_conformer(mol, cutoff_kcalmol=5.0)
        if enforce_stereo:
            smis = []
            atoms = [a.GetSymbol() for a in mol.GetAtoms()]
            for c in mol.GetConformers():
                tmp = Chem.MolFromXYZBlock(ac2xyz(atoms, c.GetPositions()))
                rdDetermineBonds.DetermineBonds(tmp)
                smis.append(Chem.MolToSmiles(tmp))
            tmp_mol = copy.deepcopy(mol)
            tmp_mol.RemoveAllConformers()
            for smi, conf in zip(smis, mol.GetConformers()):
                if smi == original_smi:
                    # only add conformers that match the original stereochemistry
                    tmp_mol.AddConformer(conf, assignId=True)
            mol = tmp_mol
        # do higher lvl single points
        sp_options = {"r2SCAN-3c": None}
        if solvent is not None:
            orca_solvent = canonicalize_solvent(solvent, "orca")
            sp_options["smd"] = orca_solvent
        r2_sp = ConformerCalculator(orca_calculate, options=sp_options, scr=scr)
        mol, _ = r2_sp(mol, n_cores=n_cores, memory=n_cores * 2)

        return mol
