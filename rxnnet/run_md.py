import argparse
import json
import logging
import os
import subprocess
import tempfile

import numpy as np
from rdkit.Chem import rdmolops, rdmolfiles
from rxnnet.mtd import md_step

from rxnnet.config import Config
from rxnnet.node import Node

DEFAULT_MD_SETTINGS = """\
$md
   temp=298.1500
   time=10
   dump=10.0000
   sdump=250
   step=0.4000
   velo=false
   shake=0
   hmass=2
   sccacc=2.0000
   nvt=true
   restart=false
$end
$metadyn
   save=250
   kpush=0.1500
   alp=0.3000
   static=false
   ramp=0.0300
$end
$wall
   potential=logfermi
   sphere:auto, all
   beta=10.0000
   temp=6000.0000
$end
$scc
   temp=6000
$end
$cma
"""


def edit_settings(settings: str) -> str:
    """Open settings in terminal editor for modification."""
    editor = os.environ.get("EDITOR", "vi")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("# MD Settings (lines starting with # are ignored)\n")
        f.write(settings)
        tmp_path = f.name

    try:
        subprocess.run([editor, tmp_path], check=True)
        with open(tmp_path, "r") as f:
            edited = f.read()
        # Remove comment lines
        lines = [
            line for line in edited.splitlines() if not line.strip().startswith("#")
        ]
        return "\n".join(lines)
    finally:
        os.unlink(tmp_path)


def run_md(
    node: Node,
    config: Config,
    multiplicity: int,
    n_cores: int,
    md_settings: str,
    run_id: str,
):
    """Run a single MD simulation for a node."""
    scr = os.getenv("SCRATCH", ".")

    # Get mol and extract atoms/coords
    lot = node.list_lots()[0] if node.list_lots() else "gfn2-xtb"
    mol = node._data[lot]
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = np.array(mol.GetConformer().GetPositions())
    charge = rdmolops.GetFormalCharge(mol)
    ranks = rdmolfiles.CanonicalRankAtoms(mol, breakTies=False)
    anum2rank = {a.GetIdx(): r for a, r in zip(mol.GetAtoms(), ranks)}

    # Build options
    solvent = config.settings.get("qm", {}).get("solvent")
    options = {"etemp": 6000}
    if solvent:
        options["alpb"] = solvent

    product = md_step(
        atoms,
        coords,
        charge,
        multiplicity,
        options=options,
        detailed_input_str=md_settings,
        n_md_cores=n_cores,
        max_products=1,
        scr=scr,
        save_traj=True,
    )
    if not product:
        raise RuntimeError("MD simulation failed to produce a product.")
    product["origin"] = [node.idx]
    origin_str = "-".join([str(idx) for idx in [node.idx]])
    product["anum2rank"] = anum2rank

    # Save product
    output_file = config.product_data / f"mtd-product_{origin_str}_{run_id}.json"
    with open(output_file, "w") as f:
        for k, v in product.items():
            if isinstance(v, np.ndarray):
                product[k] = v.tolist()
        json.dump(product, f)

    return output_file


def main():
    """Run MD simulation on a node."""
    parser = argparse.ArgumentParser(
        description="Run molecular dynamics simulation on a node"
    )
    parser.add_argument("node_file", help="Path to the node JSON file")
    parser.add_argument(
        "-d",
        "--network-dir",
        default=".",
        help="Network directory path (default: current directory)",
    )
    parser.add_argument(
        "-m", "--multiplicity", type=int, default=1, help="Multiplicity (default: 1)"
    )
    parser.add_argument(
        "-n", "--n-cores", type=int, default=1, help="Number of cores (default: 1)"
    )
    parser.add_argument(
        "-e",
        "--edit-settings",
        action="store_true",
        help="Edit MD settings in terminal editor before running",
    )
    parser.add_argument(
        "-r",
        "--run-id",
        default="0",
        help="Run identifier for output file naming (default: 0)",
    )

    args = parser.parse_args()

    config = Config(args.network_dir)

    if not config.is_initialized():
        print("Network not initialized. Run setup first.")
        return

    # Load node
    node = Node.load(args.node_file)
    print(f"Loaded node {node.idx} from {args.node_file}")

    # Get MD settings
    md_settings = DEFAULT_MD_SETTINGS
    if args.edit_settings:
        md_settings = edit_settings(md_settings)
        print("Using custom MD settings")

    # Run simulation
    print(f"Starting MD simulation with {args.n_cores} cores...")
    md_logger = logging.getLogger("rxnnet.mtd")
    md_logger.setLevel(logging.DEBUG)
    md_logger.addHandler(logging.StreamHandler())
    output_file = run_md(
        node=node,
        config=config,
        multiplicity=args.multiplicity,
        n_cores=args.n_cores,
        md_settings=md_settings,
        run_id=args.run_id,
    )
    print(f"MD simulation completed. Product saved to {output_file}")


if __name__ == "__main__":
    main()
