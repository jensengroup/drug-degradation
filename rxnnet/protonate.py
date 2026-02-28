import argparse
import copy
import os
from typing import List, Optional
from rdkit import Chem
import re
from rdkit.Chem import rdDetermineBonds, rdDistGeom, rdMolHash, rdmolops
from rxnnet.chemutils import (
    ac2mol,
    ac2xyz,
    canonicalize_solvent,
    read_multi_xyz,
    generate_and_save,
)
from rxnnet.utils import WorkingDir, check_executable, stream
from rxnnet.xtb import xtb_calculate

from rxnnet.config import Config
from rxnnet.node import Node


def get_element_graph(mol):
    """Get element graph hash for molecule comparison."""
    copy_mol = copy.deepcopy(mol)
    Chem.RemoveStereochemistry(copy_mol)
    m = Chem.RemoveAllHs(copy_mol)
    return rdMolHash.MolHash(m, rdMolHash.HashFunction.ElementGraph)


class Protonator:
    """Generate protonated or deprotonated structures using CREST."""

    def __init__(self, mode: str, crest_executable: str = "crest", swel=None):
        self.mode = mode
        if check_executable(crest_executable) is not None:
            raise RuntimeError("CREST executable not found. Please install CREST.")
        self.crest_executable = crest_executable
        self.swel = swel

    def __call__(
        self,
        mol: Chem.Mol,
        solvent: str | None = None,
        n_cores: int = 1,
    ) -> List[Chem.Mol]:
        """Run protonation/deprotonation on the given molecule."""
        scr = os.getenv("SCRATCH", ".")
        assert self.mode.lower() in [
            "protonate",
            "deprotonate",
            "add_ion",
        ], f"Invalid mode: {self.mode}. Must be 'protonate' or 'deprotonate'."

        mol = Chem.AddHs(mol)
        original_charge = rdmolops.GetFormalCharge(mol)

        def _get_charge(s: str):
            return (
                int(re.search(r"\d+", s).group()) if re.search(r"\d+", s) else 1
            ) * (-1 if s.strip().endswith("-") else 1)

        if self.mode.lower() == "protonate":
            if self.swel:
                charge = original_charge + _get_charge(self.swel)
            else:
                charge = original_charge + 1
        elif self.mode.lower() == "deprotonate":
            charge = original_charge - 1

        rdDistGeom.EmbedMolecule(mol, randomSeed=42)
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = mol.GetConformer().GetPositions()

        opt_options = {"opt": None}
        if solvent is not None:
            opt_options["alpb"] = canonicalize_solvent(solvent, qm="xtb")

        opt_results = xtb_calculate(
            atoms, coords, charge=original_charge, options=opt_options
        )
        if not opt_results["normal_termination"]:
            raise RuntimeError(f"XTB optimization failed:\n{opt_results['log']}")

        opt_coords = opt_results["opt_coords"]
        work_dir = WorkingDir(root=scr)

        with open(work_dir / "input.xyz", "w") as f:
            f.write(ac2xyz(atoms, opt_coords))

        cmd = f"{self.crest_executable} input.xyz --chrg {original_charge} --T {n_cores} --{self.mode.lower()}"
        if solvent is not None:
            cmd += f" --alpb {canonicalize_solvent(solvent, qm='xtb')}"
        if self.swel:
            cmd += f" --swel {self.swel}"
        cmd += " | tee protomers.log"

        lines = list(stream(cmd, cwd=str(work_dir)))
        normal_termination = False
        for line in reversed(lines):
            if "CREST terminated normally" in line:
                normal_termination = True
                break
        if not normal_termination:
            raise RuntimeError(f"CREST failed to run:\n{lines}")

        try:
            all_atoms, all_coords, energies = read_multi_xyz(
                work_dir / f"{self.mode.lower()}d.xyz",
                extract_property_function=lambda x: float(x),
            )
        except Exception as e:
            print(e)
            return []
        finally:
            work_dir.cleanup()

        out_mols = [ac2mol(a, c) for a, c in zip(all_atoms, all_coords)]
        for om, energy in zip(out_mols, energies):
            rdDetermineBonds.DetermineBonds(om, charge=charge)
            om.SetDoubleProp("electronic_energy", energy)

        if self.swel:
            mols = out_mols
        else:
            mols = [
                om for om in out_mols if get_element_graph(om) == get_element_graph(mol)
            ]
        return mols


def run_protonate(
    node: Node,
    config: Config,
    mode: str,
    max_charge: int,
    n_cores: int,
) -> Optional[str]:
    """Run protonation/deprotonation for a single node."""
    lot = node.list_lots()[0] if node.list_lots() else "gfn2-xtb"
    mol = node._data[lot]
    charge = rdmolops.GetFormalCharge(mol)

    # Check charge limits
    if mode == "protonate" and charge >= max_charge:
        return (
            f"Skipping node {node.idx} with charge {charge} (max allowed: {max_charge})"
        )
    if mode == "deprotonate" and charge <= -max_charge:
        return f"Skipping node {node.idx} with charge {charge} (min allowed: {-max_charge})"

    protonator = Protonator(mode=mode)
    solvent = config.settings.get("qm", {}).get("solvent")

    result = generate_and_save(
        protonator,
        {"mol": mol, "solvent": solvent, "n_cores": n_cores},
        [node.idx],
        mode.capitalize(),
        config.new_nodes,
    )
    return result


def main():
    """Generate protonated/deprotonated structures for selected nodes."""
    parser = argparse.ArgumentParser(
        description="Generate protonated or deprotonated structures for nodes"
    )
    parser.add_argument(
        "node_files",
        nargs="*",
        help="Path(s) to node JSON file(s). If not provided, interactive selection.",
    )
    parser.add_argument(
        "-d",
        "--network-dir",
        default=".",
        help="Network directory path (default: current directory)",
    )
    parser.add_argument(
        "-m",
        "--mode",
        choices=["protonate", "deprotonate"],
        default="protonate",
        help="Protonation mode (default: protonate)",
    )
    parser.add_argument(
        "-c",
        "--max-charge",
        type=int,
        default=1,
        help="Maximum absolute charge to consider (default: 1)",
    )
    parser.add_argument(
        "-n",
        "--n-cores",
        type=int,
        default=1,
        help="Number of cores to use (default: 1)",
    )

    args = parser.parse_args()

    config = Config(args.network_dir)

    if not config.is_initialized():
        print("Network not initialized. Run setup first.")
        return

    # Process nodes

    for node_file in args.node_files:
        node = Node.load(node_file)
        print(f"Processing node {node.idx} ({args.mode})...")

        _ = run_protonate(
            node=node,
            config=config,
            mode=args.mode,
            max_charge=args.max_charge,
            n_cores=args.n_cores,
        )


if __name__ == "__main__":
    main()
