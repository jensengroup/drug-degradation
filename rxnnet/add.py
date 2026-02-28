"""Add new molecules to the network from SMILES strings."""

import argparse
import os

from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdmolops
from rxnnet.xtb import xtb_calculate

from rxnnet.config import Config, IdxCounter
from rxnnet.utils import get_random_str
from rxnnet.chemutils import make_mol


def preopt(atoms, coords, charge, multiplicity, idx, solvent, n_cores, config):
    """Pre-optimization function for new molecules."""
    scr = os.getenv("SCRATCH", ".")
    # xtb optimize it
    results = xtb_calculate(
        atoms,
        coords,
        charge,
        multiplicity,
        options={"opt": None, "XTBFF": None} | ({"alpb": solvent} if solvent else {}),
        scr=scr,
        n_cores=n_cores,
    )
    results = xtb_calculate(
        atoms,
        results["opt_coords"],
        charge,
        multiplicity,
        options={"opt": None, "XTB2": None} | ({"alpb": solvent} if solvent else {}),
        scr=scr,
        n_cores=n_cores,
    )
    # make species
    mol = make_mol(results, "opt_coords")
    mol.SetProp("origin-ids", str([0]))
    mol.SetProp("origin-type", "add")
    with Chem.SDWriter(config.new_nodes / f"{get_random_str()}.sdf") as writer:
        writer.write(mol)

    return mol


def main():
    """Add new molecules to the network from SMILES strings."""
    parser = argparse.ArgumentParser(
        description="Add new molecules to the reaction network from SMILES"
    )
    parser.add_argument(
        "-d",
        "--network-dir",
        default=".",
        help="Network directory path (default: current directory)",
    )
    parser.add_argument(
        "-s",
        "--smiles",
        nargs="+",
        required=True,
        help="One or more SMILES strings for molecules",
    )
    parser.add_argument(
        "-m",
        "--multiplicity",
        type=int,
        default=1,
        help="Multiplicity of the molecule (default: 1)",
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

    node_counter = IdxCounter(str(config.network_file), "node_count")

    results = []
    for smi_str in args.smiles:
        mol = Chem.MolFromSmiles(smi_str)
        if mol is None:
            print(f"Skipping invalid SMILES: {smi_str}")
            continue
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol)
        rdDistGeom.EmbedMolecule(mol)

        # Use provided charge or auto-detect
        charge = rdmolops.GetFormalCharge(mol)
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = mol.GetConformer().GetPositions()
        node_counter.increment()

        try:
            mol = preopt(
                atoms,
                coords,
                charge,
                args.multiplicity,
                node_counter.idx,
                config.settings.get("qm", {}).get("solvent"),
                args.n_cores,
                config,
            )
            results.append(f"{smi_str}: Successfully added as node {node_counter.idx}")
        except Exception as e:
            results.append(f"{smi_str}: Failed - {e}")

    for r in results:
        print(r)


if __name__ == "__main__":
    main()
