"""Plot reactions command for visualizing reaction statistics."""

import typer
from typing_extensions import Annotated

from d2.config import NetworkConfig
from d2.utils import Reaction, choose_node

app = typer.Typer()


@app.command("plot-reactions")
def plot_reactions_command(
    node_id: Annotated[int, typer.Option(help="Index of reactant node.")] = None,
):
    """Plot reaction statistics and distributions for a given node."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    if node_id is None:
        node_files = choose_node(config.node_data, multiple=True, types=["species"])
        if node_files:
            node_id = ",".join([f.name.split("-")[0] for f in node_files])
            print(f"Selected nodes: {node_id}")
        else:
            print("No nodes selected.")
            return

    # list all reactions that start from this node
    reaction_files = list(config.reaction_data.glob(f"{node_id}-*.json"))

    if not reaction_files:
        print(f"No reactions found starting from node {node_id}")
        return

    reactions = [Reaction.load(f) for f in reaction_files]
    # sort reactions by counter
    reactions.sort(key=lambda x: x.count, reverse=True)

    # show the product distribution as a table in the terminal
    print(f"Product distribution for node {node_id}:")
    print("-" * 50)
    print(f"{'Product ID':<20}{'Reaction Type':<20}{'Count':<10}")
    print("-" * 50)
    for reaction in reactions:
        print(f"{reaction.product_id:<20}{reaction.rxn_type:<20}{reaction.count:<10}")

    print(f"\nTotal reactions: {len(reactions)}")
    print(f"Total reaction events: {sum(r.count for r in reactions)}")
