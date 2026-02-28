"""Setup command for initializing the network."""

import typer
from typing_extensions import Annotated

from d2.config import NetworkConfig

app = typer.Typer()


@app.command("setup")
def setup_command(
    temperature: Annotated[
        float,
        typer.Option(help="Temperature (optional, will prompt if not provided)"),
    ] = None,
    pressure: Annotated[
        float,
        typer.Option(help="Pressure (optional, will prompt if not provided)"),
    ] = None,
    solvent: Annotated[
        str, typer.Option(help="Solvent (optional, will prompt if not provided)")
    ] = None,
    lot: Annotated[str, typer.Option(help="Default level of theory (LOT)")] = None,
):
    """Initialize a new network with the specified parameters."""
    config = NetworkConfig()
    config.setup_network(
        temperature=temperature, pressure=pressure, solvent=solvent, lot=lot
    )
