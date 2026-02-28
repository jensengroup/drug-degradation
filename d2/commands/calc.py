"""Calc command for running calculations on nodes."""

import threading

import typer
from rdkit import Chem

from d2.calculations import optfreq
from d2.config import NetworkConfig
from d2.network import Node
from d2.utils import select_and_filter_nodes

app = typer.Typer()


@app.command("calc")
def calc_command(
    node_ids: list[str] = typer.Argument(
        None,
        help="Node ID(s) to process. Provide one or more IDs, or 'all' to select all nodes. If no IDs are provided, all nodes will be considered.",
    ),
    level: str = typer.Option("normal", help="Level of theory for the calculation"),
):
    """Run qm calculations on selected nodes."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    if node_ids and "all" in node_ids:
        node_ids = "all"
    elif node_ids:
        try:
            node_ids = [int(node_id) for node_id in node_ids]
        except ValueError:
            print(f"Invalid node IDs: {node_ids}. Must be integers or 'all'.")
            return

    node_files = select_and_filter_nodes(config, node_ids, types=["species"])
    print(f"running for {node_files}")
    if not node_files:
        return

    # skip nodes that already have a qm_data file
    filtered_files = []
    for f in node_files:
        qm_file = config.qm_data / f"{f.stem}-{level}.sdf"
        if qm_file.is_file():
            print(f"Skipping {f.stem} as it already has a QM data file.")
        else:
            filtered_files.append(f)
    node_files = filtered_files

    if not node_files:
        print("All selected nodes already have QM data.")
        return

    mols = [Node.load(f).get_mol(config.lot) for f in node_files]
    sdfs = [Chem.MolToMolBlock(mol) for mol in mols]
    names = [f.stem for f in node_files]
    multiplicity = 1

    print(f"Starting calculations for {len(node_files)} nodes...")

    for name, sdf in zip(names, sdfs):
        t = threading.Thread(
            target=optfreq,
            args=(
                sdf,
                name,
                level,
                multiplicity,
                config.solvent,
                1.0,
                str(config.qm_data),
            ),
            daemon=False,
        )
        t.start()
        print(f"Started calculation for {name}")
