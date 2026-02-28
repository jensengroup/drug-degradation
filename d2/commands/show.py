"""Show command for viewing molecular structures."""

import webbrowser
from pathlib import Path

import typer
from rdkit import Chem
from rich.columns import Columns
from rich.panel import Panel
from rich.prompt import Prompt
from typing_extensions import Annotated

from d2.config import NetworkConfig
from d2.network import Node

app = typer.Typer()


def choose_node(dir: Path, multiple: bool = False, types=["species", "complex"]):
    """Choose node files interactively."""
    files = [f for f in dir.glob("*.json")]
    if not files:
        print("No nodes available in the network.")
        return
    species_files = sorted([f.name for f in files if int(f.name.split("-")[0]) < 1000])
    complex_files = sorted([f.name for f in files if int(f.name.split("-")[0]) > 999])

    # Create species and complex panels with files listed below them
    species_panel = Panel("\\n".join(species_files), title="Species")
    complex_panel = Panel("\\n".join(complex_files), title="Complexes")

    type2column = {"species": species_panel, "complex": complex_panel}

    # Create columns to display both panels side by side
    columns = Columns([type2column[t] for t in types])

    print(columns)
    if multiple:
        prompt = (
            "Select one or multiple species by index, separated by space (e.g. '1 2 3')"
        )
    else:
        prompt = "Select a species by index"
    user_input = Prompt.ask(prompt)
    if multiple:
        try:
            selected_ids = [int(idx) for idx in user_input.split(" ")]
        except ValueError:
            print(f"Invalid index provided: {user_input}")
            return
        filenames = []
        for idx in selected_ids:
            filename = [f for f in files if int(f.name.split("-")[0]) == idx]
            if len(filename) == 1:
                filenames.append(filename[0])
            else:
                print(f"Invalid index provided: {idx}")
        return filenames
    else:
        try:
            selected_id = int(user_input)
        except ValueError:
            print(f"Invalid index provided: {user_input}")
            return
        filename = [f for f in files if int(f.name.split("-")[0]) == selected_id]
        if len(filename) == 1:
            return filename[0]
        else:
            print(f"Invalid index provided: {user_input}")


def show_molecule(filename, config: NetworkConfig):
    """Load a JSON file from node_data and display its molecular structure
    using 3Dmol.js."""
    file_path = config.node_data / filename

    if not file_path.is_file():
        print(f"Error: {filename} does not exist in {config.node_data}")
        return

    n = Node.load(file_path)
    sdf = Chem.MolToMolBlock(n.get_mol(n.list_lots()[0]))

    html_content = f"""
    <html>
    <head>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    </head>
    <body>
        <div id="viewer" style="width: 100%; height: 100%;"></div>
        <script>
            let viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "white" }});
            let sdfData = `{sdf}`;
            viewer.addModel(sdfData, "sdf");
            viewer.setStyle({{"sphere": {{"radius": 0.4}}, "stick": {{}}}})
            viewer.setBackgroundColor("0xeeeeee", 0)

            viewer.zoomTo();
            viewer.render();
        </script>
    </body>
    </html>
    """

    html_file = config.network_dir / ".molecule_viewer.html"
    with open(html_file, "w") as f:
        f.write(html_content)

    webbrowser.open(f"file://{html_file}")


@app.command("show")
def show_command(
    filename: Annotated[
        str, typer.Option(help="Filename (optional, will prompt if not provided)")
    ] = None,
):
    """Select a node file from node_data and view it using 3Dmol.js."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    if filename is None:
        filename = choose_node(config.node_data)
    if filename:
        show_molecule(filename, config)
