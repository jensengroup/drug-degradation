"""Protonate command for generating protonated structures."""

import argparse
from typing import Optional
from rxnnet.chemutils import (
    generate_and_save,
)
from rxnnet.protonate import Protonator

from rxnnet.config import Config
from rxnnet.node import Node


def add_ion(
    node: Node,
    config: Config,
    ion: str,
    n_cores: int,
) -> Optional[str]:
    """Run protonation/deprotonation for a single node."""
    lot = node.list_lots()[0] if node.list_lots() else "gfn2-xtb"
    mol = node._data[lot]

    AddIonGenerator = Protonator(mode="protonate", swel=ion)
    solvent = config.settings.get("qm", {}).get("solvent")

    result = generate_and_save(
        AddIonGenerator,
        {"mol": mol, "solvent": solvent, "n_cores": n_cores},
        [node.idx],
        f"add-{ion.capitalize()}",
        config.new_nodes,
    )
    return result


def main():
    """Generate structures for selected nodes with Ion added."""
    parser = argparse.ArgumentParser(
        description="Generate structures for selected nodes with Ion added"
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
        "-i",
        "--ion",
        type=str,
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
        print(f"Processing node {node.idx}...")

        out = add_ion(
            node=node,
            config=config,
            ion=args.ion,
            n_cores=args.n_cores,
        )
        print(out)


if __name__ == "__main__":
    main()
