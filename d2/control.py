import hashlib
import json
import logging
import math
import os
import re
import threading
import webbrowser
from itertools import product
from pathlib import Path

import click
import numpy as np
import submitit
import typer
from pyvis.network import Network
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor, rdDetermineBonds, rdDistGeom, rdmolops
from rich import print
from rich.columns import Columns
from rich.panel import Panel
from rich.prompt import Prompt
from tooltoad.chemutils import ac2mol, hartree2kcalmol
from tooltoad.mdd import md_step
from tooltoad.orca import orca_calculate
from tooltoad.xtb import (
    MDOptions,
    MetaDynOptions,
    SCCOptions,
    WallOptions,
    xtb_calculate,
)
from typing_extensions import Annotated

from d2.calculations import optfreq
from d2.network import Node
from d2.protonation import Protonator
from d2.tautomers import TautomerGenerator

OFFSETS = {
    "water": {"deprotonation": -0.4059596304589945, "protonation": 0.4059596304589945}
}

app = typer.Typer()

network_dir = "."

SCR = os.getenv("SCRATCH", ".")

network_dir = Path(network_dir).resolve()
network_file = network_dir / "network-info.json"
node_data = network_dir / "nodes"
product_data = network_dir / "products"
new_nodes = network_dir / "new-nodes"
reaction_data = network_dir / "reactions"
qm_data = network_dir / "qm-data"

# load global settings
if network_file.is_file():
    with open(network_file, "r") as f:
        global_settings = json.load(f)
    TEMPERATURE = global_settings["temperature"]
    PRESSURE = global_settings["pressure"]
    SOLVENT = global_settings.get("solvent", None)
    NODE_COUNT = global_settings.get("node_count", 0)
    LOT = global_settings.get("lot", "gfn2-xtb")


else:
    print("No `network-info.json` file found, run `control.py setup`")
    NODE_COUNT = 0
#    sys.exit()
# Set up logging


executor = submitit.AutoExecutor(
    folder=".tmp",
)
executor.update_parameters(
    name="md",
    cpus_per_task=4,
    timeout_min=1200,
    slurm_partition="kemi1",
    slurm_mem_per_cpu=4000,
    slurm_array_parallelism=100,
)


class IdxCounter:
    def __init__(self, parent_file: str, name: str):
        self.parent_file = parent_file
        self.name = name

    @property
    def idx(self):
        with open(self.parent_file, "r") as f:
            global_settings = json.load(f)
        return global_settings.get(self.name)

    def increment(self):
        with open(self.parent_file, "r") as f:
            global_settings = json.load(f)
        global_settings[self.name] += 1
        with open(self.parent_file, "w") as f:
            json.dump(global_settings, f, indent=2)
        return global_settings[self.name]


node_counter = IdxCounter(network_file, "node_count")


@app.command()
def setup(
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
    solvent_choices = [
        "None",
        "Acetone",
        "Acetonitrile",
        "CH2Cl2",
        "CHCl3",
        "DMF",
        "DMSO",
        "Hexane",
        "Methanol",
        "Octanol",
        "THF",
        "Water",
    ]
    if solvent is not None and solvent not in solvent_choices:
        raise typer.BadParameter(
            f"Invalid solvent choice. Must be one of {solvent_choices}"
        )
    # Check if the network file exists
    if network_file.is_file():
        print(f"Network is already initialized in {network_dir}.")
        return

    print(f"No Network initialized in {network_dir}.")
    if not typer.confirm("Do you want to initialize a new Network?", default=True):
        return

    if lot is None:
        lot = Prompt.ask(
            "Enter default level of theory (LOT)",
            default="r2scan-3c/gfn2-xtb",
            show_default=True,
        ).lower()

    # Prompt only for missing values
    if temperature is None:
        temperature = Prompt.ask("Enter temperature (K)", default=298.15)
    if pressure is None:
        pressure = Prompt.ask("Enter pressure (atm)", default=1.0)
    if solvent is None:
        solvent = Prompt.ask(
            "Enter solvent",
            choices=solvent_choices,
            default="None",
        )

    # Proceed with setup
    print(
        f"Initializing network with Temperature={temperature}, Pressure={pressure}, Solvent={solvent}"
    )

    # Store global settings in a JSON file
    global_settings = (
        {
            "temperature": temperature,
            "pressure": pressure,
            "node_count": 0,
            "lot": lot,
        }
        | {"solvent": solvent.lower()}
        if solvent != "None"
        else {}
    )

    with open(network_file, "w") as f:
        json.dump(global_settings, f, indent=2)
    os.makedirs(node_data.resolve())
    os.makedirs(reaction_data.resolve())
    os.makedirs(product_data.resolve())
    os.makedirs(new_nodes.resolve())
    os.makedirs(qm_data.resolve())
    print(f"Network initialized at {network_file}")

    # be able to select a node and view it
    # a method called show which takes a filename of a json file in node_data and loads its contecnt and shows a molecule uing 3dmol.js


# @app.command()
# def vis(
#     substrate_id: int = 1, lot: str = "gfn2-xtb", prop_name: str = "electronic_energy"
# ):
#     ENERGY_TYPE = "G" if "/" in lot else "E"
#     rdDepictor.SetPreferCoordGen(True)
#     nodes = [Node.load(f) for f in node_data.glob("*.json")]
#     # make maps
#     mol_svg_map = {}
#     mol_energy_map = {}
#     for n in nodes:
#         mol = n.get_mol(lot)
#         mol_energy_map[n.idx] = mol.GetDoubleProp(prop_name)

#         try:
#             mol = Chem.RemoveHs(mol)
#         except Exception as e:
#             pass
#         mol.RemoveAllConformers()
#         img = Draw.MolsToGridImage(
#             [
#                 mol,
#             ],
#             molsPerRow=1,
#             subImgSize=(500, 400),
#             useSVG=True,
#         )

#         mol_svg_map[n.idx] = img

#     # now get all the reactions
#     reactions = [Reaction.load(f) for f in reaction_data.glob("*.json")]

#     edges = [
#         {
#             "begin": r.reactant_id,
#             "end": r.product_id,
#             "type": r.rxn_type,
#             "deltaE": hartree2kcalmol(
#                 mol_energy_map.get(r.product_id, 0)
#                 - mol_energy_map.get(r.reactant_id, 0)
#                 + OFFSETS.get(SOLVENT.lower(), {}).get(r.rxn_type.lower(), 0)
#             ),
#         }
#         for r in reactions
#     ]
#     # === Build the PyVis network ===
#     net = Network(
#         height="800px",
#         width="1200px",
#         directed=True,
#         notebook=True,
#         cdn_resources="in_line",
#     )

#     # === Add nodes ===
#     for n in nodes:
#         if n.idx == substrate_id:
#             net.add_node(
#                 n.idx,
#                 label="Substrate",
#                 shape="ellipse",  # fixed typo: "elipse" → "ellipse"
#                 color="#4E4E4EFF",
#                 font={"color": "white"},
#             )
#         else:
#             net.add_node(
#                 n.idx,
#                 label=str(n.idx),
#                 shape="circle",
#                 color="#C2C2C2",
#             )

#     # === Add edges ===
#     for edge in edges:
#         net.add_edge(
#             edge["begin"],
#             edge["end"],
#             title=f"{edge['type']}\nΔ{ENERGY_TYPE} {edge['deltaE']:.02f} kcal/mol",
#             color="#8B29E184",  # consider edge color per type here if needed
#         )

#     # === Export to base HTML ===
#     base_filename = ".reaction_network_vis.html"
#     net.show(base_filename)

#     # === Custom HTML UI: Slider + Modal Image Preview ===
#     modal_html = f"""
#     <!-- UI Styling -->
#     <style>
#     html, body {{
#         margin: 0;
#         padding: 0;
#         height: 100%;
#         width: 100%;
#         background: white;
#         overflow: hidden;
#     }}
#     #mynetwork {{
#         width: 100%;
#         height: 100%;
#         border: none;
#         background: none;
#         box-shadow: none;
#     }}
#     canvas {{
#         background: transparent;
#         border: none;
#         display: block;
#         width: 100%;
#         height: 100%;
#     }}
#     .vis-network {{
#         width: 100%;
#         height: 100%;
#     }}
#     </style>

#     <!-- Threshold slider -->
#     <div style="position: fixed; top: 10px; left: 10px; background: white; padding: 10px; border: 1px solid #ccc; z-index: 9999;">
#         <label for="threshold">Edge Count Threshold: <span id="thresholdValue">5</span></label>
#         <input type="range" id="threshold" min="1" max="10" value="5" style="width: 200px;">
#     </div>

#     <!-- Modal for displaying molecule SVG -->
#     <div id="imageModal" style="
#         position: fixed;
#         top: 80px;
#         left: 20px;
#         background: white;
#         border: 2px solid #444;
#         padding: 15px;
#         display: none;
#         z-index: 9999;
#         box-shadow: 0px 0px 10px #aaa;
#     ">
#         <div id="modalImage" style="
#             width: 100%;
#             max-width: 500px;
#             margin-bottom: 10px;
#             display: flex;
#             align-items: center;
#             justify-content: center;
#         "></div>
#         <button onclick="document.getElementById('imageModal').style.display='none';">Close</button>
#     </div>

#     <!-- Interactive logic -->
#     <script type="text/javascript">
#         const svgMap = {json.dumps(mol_svg_map)};

#         // Threshold slider
#         document.getElementById("threshold").addEventListener("input", function () {{
#             const val = parseInt(this.value);
#             document.getElementById("thresholdValue").innerText = val;
#             updateEdges(val);
#         }});

#         // Image popup on node click
#         network.on("click", function (params) {{
#             if (params.nodes.length > 0) {{
#                 const nodeId = params.nodes[0];
#                 const svg = svgMap[nodeId];
#                 if (svg) {{
#                     document.getElementById("modalImage").innerHTML = svg;
#                     document.getElementById("imageModal").style.display = "block";
#                 }} else {{
#                     console.warn("No SVG found for node:", nodeId);
#                 }}
#             }}
#         }});

#         // Initialize edges based on current threshold
#         updateEdges(parseInt(document.getElementById("threshold").value));
#     </script>
#     """

#     # === Inject modal HTML before </body> in the PyVis-generated HTML ===
#     html_path = Path(base_filename)
#     html_content = html_path.read_text(encoding="utf-8")
#     modified_html = html_content.replace("</body>", modal_html + "\n</body>")
#     html_path.write_text(modified_html, encoding="utf-8")

#     # open the file in the browser
#     webbrowser.open(f"file://{html_path.resolve()}")


@app.command()
def vis(
    substrate_id: int = 1, lot: str = "gfn2-xtb", prop_name: str = "l2l1_gibbs-energy"
):
    # ENERGY_TYPE = "G" if "/" in lot else "E"
    ENERGY_TYPE = "G" if "gibbs" in prop_name else "E"
    rdDepictor.SetPreferCoordGen(True)
    nodes = [Node.load(f) for f in node_data.glob("*.json")]

    mol_svg_map = {}
    mol_energy_map = {}
    for n in nodes:
        # mol = n.get_mol(lot)
        mol = n.get_qm(qm_data)
        if mol is None:
            print(f"Node {n.idx} has no QM data, skipping.")
            continue
        mol_energy_map[n.idx] = mol.GetDoubleProp(prop_name)

        try:
            mol = Chem.RemoveHs(mol)
        except Exception:
            pass
        mol.RemoveAllConformers()
        img = Draw.MolsToGridImage(
            [mol], molsPerRow=1, subImgSize=(500, 400), useSVG=True
        )
        mol_svg_map[n.idx] = img

    substrate_energy = mol_energy_map[substrate_id]
    rel_energy_map = {
        idx: energy - substrate_energy for idx, energy in mol_energy_map.items()
    }

    reactions = [Reaction.load(f) for f in reaction_data.glob("*.json")]
    edges = [
        {
            "begin": r.reactant_id,
            "end": r.product_id,
            "type": r.rxn_type,
            "deltaE": hartree2kcalmol(
                mol_energy_map.get(r.product_id, 0)
                - mol_energy_map.get(r.reactant_id, 0)
                + OFFSETS.get(SOLVENT.lower(), {}).get(r.rxn_type.lower(), 0)
            ),
        }
        for r in reactions
    ]

    net = Network(
        height="800px",
        width="1200px",
        directed=True,
        notebook=True,
        cdn_resources="in_line",
    )

    for n in nodes:
        color = "#4E4E4EFF" if n.idx == substrate_id else "#C2C2C2"
        shape = "ellipse" if n.idx == substrate_id else "circle"
        label = "Substrate" if n.idx == substrate_id else str(n.idx)
        font = {"color": "white"} if n.idx == substrate_id else {}

        net.add_node(n.idx, label=label, shape=shape, color=color, font=font)

    for i, edge in enumerate(edges):
        net.add_edge(
            edge["begin"],
            edge["end"],
            title=f"{edge['type']}\n∆{ENERGY_TYPE} {edge['deltaE']:.02f} kcal/mol",
            color="#8B29E184",
            id=i,
        )

    base_filename = ".reaction_network_vis.html"
    net.show(base_filename)

    modal_html = f"""
    <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>
    <style>
      html, body {{ margin: 0; padding: 0; height: 100%; width: 100%; background: white; overflow: hidden; }}
      #mynetwork, .vis-network, canvas {{ width: 100%; height: 100%; border: none; background: none; box-shadow: none; }}
    </style>

    <div style=\"position: fixed; top: 10px; left: 10px; background: white; padding: 10px; border: 1px solid #ccc; z-index: 9999;\">
        <label for=\"threshold\">Edge Count Threshold: <span id=\"thresholdValue\">5</span></label>
        <input type=\"range\" id=\"threshold\" min=\"1\" max=\"10\" value=\"5\" style=\"width: 200px;\">
    </div>

    <div style=\"position: fixed; bottom: 10px; left: 10px; background: white; padding: 10px; border: 1px solid #ccc; z-index: 9999;\">
        <div id=\"histogram\" style=\"width: 300px; height: 200px;\"></div>
        <label for=\"energyCutoff\">Energy Cutoff: <span id=\"energyCutoffValue\">20</span> kcal/mol</label>
        <input type=\"range\" id=\"energyCutoff\" min=\"0\" max=\"100\" step=\"1\" value=\"20\" style=\"width: 300px;\">
    </div>

    <div id=\"imageModal\" style=\"position: fixed; top: 80px; left: 20px; background: white; border: 2px solid #444; padding: 15px; display: none; z-index: 9999; box-shadow: 0px 0px 10px #aaa;\">
        <div id=\"modalImage\" style=\"width: 100%; max-width: 500px; margin-bottom: 10px; display: flex; align-items: center; justify-content: center;\"></div>
        <button onclick=\"document.getElementById('imageModal').style.display='none';\">Close</button>
    </div>

    <script type=\"text/javascript\">
        const svgMap = {json.dumps(mol_svg_map)};
        const relEnergyMap = {json.dumps(rel_energy_map)};
        let visibleNodes = new Set(Object.keys(relEnergyMap).map(k => parseInt(k)));

        function updateEnergyFilter(cutoff) {{
            const nodes = network.body.data.nodes;
            const edges = network.body.data.edges;
            visibleNodes.clear();
            nodes.forEach(n => {{
                const energy = relEnergyMap[n.id];
                const visible = energy <= cutoff;
                visibleNodes.add(n.id);
                nodes.update({{id: n.id, hidden: !visible}});
            }});
            edges.forEach(e => {{
                const visible = visibleNodes.has(e.from) && visibleNodes.has(e.to);
                edges.update({{id: e.id, hidden: !visible}});
            }});
        }}

        Plotly.newPlot("histogram", [{{
            x: Object.values(relEnergyMap),
            type: "histogram",
            marker: {{color: "#8888ff"}},
        }}], {{
            margin: {{t: 10, b: 30, l: 30, r: 10}},
            xaxis: {{title: "ΔE relative to substrate"}},
            yaxis: {{title: "Count"}},
        }});

        document.getElementById("energyCutoff").addEventListener("input", function () {{
            const val = parseFloat(this.value);
            document.getElementById("energyCutoffValue").innerText = val;
            updateEnergyFilter(val);
        }});

        network.on("click", function (params) {{
            if (params.nodes.length > 0) {{
                const nodeId = params.nodes[0];
                const svg = svgMap[nodeId];
                if (svg) {{
                    document.getElementById("modalImage").innerHTML = svg;
                    document.getElementById("imageModal").style.display = "block";
                }} else {{
                    console.warn("No SVG found for node:", nodeId);
                }}
            }}
        }});
    </script>
    """

    html_path = Path(base_filename)
    html_content = html_path.read_text(encoding="utf-8")
    modified_html = html_content.replace("</body>", modal_html + "\n</body>")
    html_path.write_text(modified_html, encoding="utf-8")
    webbrowser.open(f"file://{html_path.resolve()}")


@app.command()
def collect():
    for f in new_nodes.glob("*.sdf"):
        # check if hash already exists in the node_data directory
        mol_hash = f.stem.split("-")[0]
        suppl = Chem.SDMolSupplier(str(f), removeHs=False)
        mol = next(suppl, None)
        product_node_files = list(node_data.glob(f"*-{mol_hash}.json"))
        if len(product_node_files) > 0:
            print(f"Node {mol_hash} already exists")
            node = Node.load(product_node_files[0])
        else:
            print("need to make new node")
            node_counter.increment()
            node = Node(
                idx=node_counter.idx,
                data={"gfn2-xtb": mol},
            )
            node.save(node_data)
        origin_idx = mol.GetIntProp("origin-idx")
        origin_type = mol.GetProp("origin-type")
        # check for self reaction
        if origin_idx == node.idx:
            print(f"Skipping self-reaction for node {node.idx}")
            f.unlink()
            continue
        # check if reaction connecting origin-idx and new node already exists
        reaction_files = list(reaction_data.glob(f"{origin_idx}-{node.idx}-*.json"))
        if len(reaction_files) == 0:
            print(f"Creating reaction from {origin_idx} to {node.idx}")
            rxn = Reaction(
                reactant_id=origin_idx, product_id=node.idx, rxn_type=origin_type
            )
            rxn.save(reaction_data)
        elif origin_type == "md-reaction":
            raise NotImplementedError("needs to increment counter")
        else:
            print(
                f"nothing to do, rxn of type {origin_type} between {origin_idx} and {node.idx} already exists"
            )
        f.unlink()  # remove the file after processing


@app.command()
def collect_qm():
    for f in qm_data.glob("*.sdf"):
        pass
        # node = Node.load(f)


@app.command()
def calc(
    node_ids: Annotated[
        str,
        typer.Argument(
            help="Node ID(s) to process. Provide a single ID, space- or comma-separated list, or 'all'."
        ),
    ] = None,
    level: Annotated[
        str, typer.Option(help="Level of theory for the calculation")
    ] = "normal",
):
    if node_ids is None:
        node_files = [
            Path(f) for f in choose_node(node_data, multiple=True, types=["species"])
        ]
    elif node_ids.lower() == "all":
        node_files = list(node_data.glob("*.json"))
    else:
        raise NotImplementedError()
        node_files = [
            Node.load((node_data.glob(f"{idx}-").__next__())) for idx in node_ids
        ]

    if not node_files:
        print("No valid node files found.")
        return

    # skip nodes that already have a qm_data file
    filtered_files = []
    for f in node_files:
        qm_file = qm_data / f"{f.stem}-{level}.sdf"
        if qm_file.is_file():
            print(f"Skipping {f.stem} as it already has a QM data file.")
        else:
            filtered_files.append(f)
    node_files = filtered_files
    mols = [Node.load(f).get_mol(LOT) for f in node_files]
    sdfs = [Chem.MolToMolBlock(mol) for mol in mols]
    names = [f.stem for f in node_files]
    multiplicity = 1

    for name, sdf in zip(names, sdfs):
        t = threading.Thread(
            target=optfreq,
            args=(sdf, name, level, multiplicity, SOLVENT, 1.0, str(qm_data)),
            daemon=False,  # optional: makes threads exit when main program exits
        )
        t.start()


@app.command()
def tautomerize(
    node_ids: Annotated[
        str,
        typer.Argument(
            help="Node ID(s) to process. Provide a single ID, space- or comma-separated list, or 'all'."
        ),
    ] = None,
    n_cores: Annotated[
        int, typer.Option(help="Number of cores to use for the calculation")
    ] = 1,
    remote: Annotated[
        bool, typer.Option(help="Submit QM calculation to remote executer")
    ] = False,
):
    if node_ids is None:
        node_files = [
            Path(f) for f in choose_node(node_data, multiple=True, types=["species"])
        ]
    elif node_ids.lower() == "all":
        node_files = list(node_data.glob("*.json"))
    else:
        raise NotImplementedError()
        node_files = [
            Node.load((node_data.glob(f"{idx}-").__next__())) for idx in node_ids
        ]

    if not node_files:
        print("No valid node files found.")
        return

    tautomer_generator = TautomerGenerator()
    if remote:
        executor.update_parameters(name="tautomers", cpus_per_task=n_cores)
        jobs = []
        with executor.batch():
            for file in node_files:
                click.echo(f"Generating tautomers for: {file.stem}")
                node = Node.load(file)
                mol = node.get_mol(LOT)
                job = executor.submit(
                    generate_and_save,
                    tautomer_generator,
                    {"mol": mol, "solvent": SOLVENT, "n_cores": n_cores},
                    node.idx,
                    "Tautomer",
                )
                jobs.append(job)
        print(f"Submitted {len(jobs)} jobs for tautomers generation.")

    else:
        for file in node_files:
            click.echo(f"Generating tautomers for: {file.stem}")
            node = Node.load(file)
            mol = node.get_mol(LOT)
            result = generate_and_save(
                tautomer_generator,
                {"mol": mol, "solvent": SOLVENT, "n_cores": n_cores},
                node.idx,
                "Tautomer",
            )
            print(result)


@app.command()
def protonate(
    node_ids: Annotated[
        str,
        typer.Argument(
            help="Node ID(s) to process. Provide a single ID, space- or comma-separated list, or 'all'."
        ),
    ] = None,
    n_cores: Annotated[
        int, typer.Option(help="Number of cores to use for the calculation")
    ] = 1,
    remote: Annotated[
        bool, typer.Option(help="Submit QM calculation to remote executer")
    ] = False,
    max_charge: Annotated[int, typer.Option(help="Maximum charge to consider")] = 1,
):
    if node_ids is None:
        node_files = [
            Path(f) for f in choose_node(node_data, multiple=True, types=["species"])
        ]
    elif node_ids.lower() == "all":
        node_files = list(node_data.glob("*.json"))
    else:
        raise NotImplementedError()
        node_files = [
            Node.load((node_data.glob(f"{idx}-").__next__())) for idx in node_ids
        ]

    if not node_files:
        print("No valid node files found.")
        return

    protonator = Protonator(mode="protonate")
    if remote:
        executor.update_parameters(name="protonation", cpus_per_task=n_cores)
        jobs = []
        with executor.batch():
            for file in node_files:
                click.echo(f"Generating Protonated Structure for: {file.stem}")
                node = Node.load(file)
                mol = node.get_mol(LOT)
                charge = rdmolops.GetFormalCharge(mol)
                if charge >= max_charge:
                    print(
                        f"Skipping {file.stem} with charge {charge} (max allowed: {max_charge})"
                    )
                    continue
                job = executor.submit(
                    generate_and_save,
                    protonator,
                    {"mol": mol, "solvent": SOLVENT, "n_cores": n_cores},
                    node.idx,
                    "Protonation",
                )
                jobs.append(job)
        print(f"Submitted {len(jobs)} jobs for Protonation generation.")

    else:
        for file in node_files:
            click.echo(f"Generating Protonated Structure for: {file.stem}")
            node = Node.load(file)
            mol = node.get_mol(LOT)
            charge = rdmolops.GetFormalCharge(mol)
            if charge >= max_charge:
                print(
                    f"Skipping {file.stem} with charge {charge} (max allowed: {max_charge})"
                )
                continue
            result = generate_and_save(
                protonator,
                {"mol": mol, "solvent": SOLVENT, "n_cores": n_cores},
                node.idx,
                "Protonation",
            )
            print(result)


@app.command()
def deprotonate(
    node_ids: Annotated[
        str,
        typer.Argument(
            help="Node ID(s) to process. Provide a single ID, space- or comma-separated list, or 'all'."
        ),
    ] = None,
    n_cores: Annotated[
        int, typer.Option(help="Number of cores to use for the calculation")
    ] = 1,
    remote: Annotated[
        bool, typer.Option(help="Submit QM calculation to remote executer")
    ] = False,
    min_charge: Annotated[int, typer.Option(help="Maximum charge to consider")] = -1,
):
    if node_ids is None:
        node_files = [
            Path(f) for f in choose_node(node_data, multiple=True, types=["species"])
        ]
    elif node_ids.lower() == "all":
        node_files = list(node_data.glob("*.json"))
    else:
        raise NotImplementedError()
        node_files = [
            Node.load((node_data.glob(f"{idx}-").__next__())) for idx in node_ids
        ]

    if not node_files:
        print("No valid node files found.")
        return

    deprotonator = Protonator(mode="deprotonate")
    if remote:
        executor.update_parameters(name="deprotonation", cpus_per_task=n_cores)
        jobs = []
        with executor.batch():
            for file in node_files:
                click.echo(f"Generating Deprotonated Structure for: {file.stem}")
                node = Node.load(file)
                mol = node.get_mol(LOT)
                charge = rdmolops.GetFormalCharge(mol)
                if charge <= min_charge:
                    print(
                        f"Skipping {file.stem} with charge {charge} (max allowed: {min_charge})"
                    )
                    continue
                job = executor.submit(
                    generate_and_save,
                    deprotonator,
                    {"mol": mol, "solvent": SOLVENT, "n_cores": n_cores},
                    node.idx,
                    "Deprotonation",
                )
                jobs.append(job)
        print(f"Submitted {len(jobs)} jobs for Deprotonation generation.")

    else:
        for file in node_files:
            click.echo(f"Generating Deprotonated Structure for: {file.stem}")
            node = Node.load(file)
            mol = node.get_mol(LOT)
            charge = rdmolops.GetFormalCharge(mol)
            if charge <= min_charge:
                print(
                    f"Skipping {file.stem} with charge {charge} (max allowed: {min_charge})"
                )
                continue
            result = generate_and_save(
                deprotonator,
                {"mol": mol, "solvent": SOLVENT, "n_cores": n_cores},
                node.idx,
                "Deprotonation",
            )
            print(result)


# have command to make new dft calcs for nodes
# also one to make conf search for nodes
# to pipeline would be
# add node
# generate new nodes
# etcetc
# for energetics
# conf search each node
# dft each node


def get_random_str():
    """Generate a random string of 6 characters."""
    return hashlib.md5(os.urandom(16)).hexdigest()[:6]


def generate_and_save(generate_function, kwargs, origin: int, origin_type: str):
    results = generate_function(**kwargs)
    for result in results:
        mol_hash = hashlib.md5(Chem.MolToSmiles(result).encode()).hexdigest()[:12]
        result.SetIntProp("origin-idx", origin)
        result.SetProp("origin-type", origin_type)
        with Chem.SDWriter(new_nodes / f"{mol_hash}-{get_random_str()}.sdf") as writer:
            writer.write(result)
    return f"Generated {len(results)} molecules from {origin}."


@app.command()
def plot_reactions(
    node_id: Annotated[int, typer.Option(help="Index of reactant node.")] = None,
):
    if node_id is None:
        node_files = choose_node(node_data, multiple=True, types=["species"])
        node_id = ",".join([f.split("-")[0] for f in node_files])
        print(node_id)
    # list all reactions that start from this node
    reaction_files = list(reaction_data.glob(f"{node_id}-*.json"))
    # load all reactions
    reactions = [Reaction.load(f) for f in reaction_files]
    # sort reactions by counter
    reactions.sort(key=lambda x: x.count, reverse=True)
    # show the product distribution as a table in the terminal
    print(f"Product distribution for node {node_id}:")
    print("-" * 50)
    print(f"{'Reaction':<20}{'Count':<10}")
    print("-" * 50)
    for reaction in reactions:
        print(f"{reaction.product_id:<20}{reaction.count:<10}")


@app.command()
def show(
    filename: Annotated[
        str, typer.Option(help="Filename (optional, will prompt if not provided)")
    ] = None,
):
    """Select a node file from node_data and view it using 3Dmol.js."""
    if filename is None:
        filename = choose_node(node_data)
    if filename:
        show_molecule(filename)


def choose_node(dir: Path, multiple: bool = False, types=["species", "complex"]):
    files = [f for f in dir.glob("*.json")]
    if not files:
        print("No nodes available in the network.")
        return
    species_files = sorted([f.name for f in files if int(f.name.split("-")[0]) < 1000])
    complex_files = sorted([f.name for f in files if int(f.name.split("-")[0]) > 999])

    # Create species and complex panels with files listed below them
    species_panel = Panel("\n".join(species_files), title="Species")
    complex_panel = Panel("\n".join(complex_files), title="Complexes")

    type2column = {"species": species_panel, "complex": complex_panel}

    # Create columns to display both panels side by side
    columns = Columns([type2column[t] for t in types])

    print(columns)
    if multiple:
        prompt = (
            "Select one or multiple species by index, seperated by space (e.g. '1 2 3')"
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


def show_molecule(filename):
    """Load a JSON file from node_data and display its molecular structure
    using 3Dmol.js."""
    file_path = node_data / filename

    if not file_path.is_file():
        print(f"Error: {filename} does not exist in {node_data}")
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

    html_file = network_dir / ".molecule_viewer.html"
    with open(html_file, "w") as f:
        f.write(html_content)

    webbrowser.open(f"file://{html_file}")


# # species must have atoms and coords attribute, coords always sorted by lowest energy, the atom number can be whatever, so need to store as pairs of atoms and coords
# # have property that will always give the lowest energy atoms and coords out


def make_mol(data: dict, coords: str = "opt_coords"):
    mol = ac2mol(data["atoms"], data[coords])
    rdDetermineBonds.DetermineBonds(mol, charge=data["charge"])
    for prop in ["electronic_energy", "gibbs_energy"]:
        if prop in data:
            mol.SetDoubleProp(prop, data[prop])
    return mol


@app.command()
def add(
    smi: Annotated[
        list[str], typer.Argument(help="One or more SMILES strings for molecules")
    ],
    multiplicity: Annotated[int, typer.Option(help="Multiplicity of the Molecule")] = 1,
    remote: Annotated[
        bool, typer.Option(help="Submit QM calculation to remote executer")
    ] = False,
    n_cores: Annotated[int, typer.Option(help="Number of cores to use")] = 1,
):
    # when i use the once from jobs i need to make sure i have the options correctly set
    def preopt(atoms, coords, charge, multiplicity, idx, SOLVENT, n_cores):
        SCR = os.getenv("SCRATCH", ".")
        # xtb optimize it
        results = orca_calculate(
            atoms,
            coords,
            charge,
            multiplicity,
            options={"opt": None, "XTBFF": None}
            | ({"alpb": SOLVENT} if SOLVENT else {}),
            scr=SCR,
            n_cores=n_cores,
        )
        results = orca_calculate(
            atoms,
            results["opt_coords"],
            charge,
            multiplicity,
            options={"opt": None, "XTB2": None}
            | ({"alpb": SOLVENT} if SOLVENT else {}),
            scr=SCR,
            n_cores=n_cores,
        )
        # make species
        mol = make_mol(results, "opt_coords")
        node = Node(
            data={"gfn2-xtb": mol},
            idx=idx,
        )
        node.save(node_data)
        return "Done"

    results = []
    for smi_str in smi:
        mol = Chem.MolFromSmiles(smi_str)
        if mol is None:
            print(f"Skipping invalid SMILES: {smi_str}")
            continue
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol)
        rdDistGeom.EmbedMolecule(mol)
        charge = rdmolops.GetFormalCharge(mol)
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = mol.GetConformer().GetPositions()
        node_counter.increment()

        if remote:
            executor.update_parameters(name="preopt-species", cpus_per_task=n_cores)
            job = executor.submit(
                preopt,
                atoms,
                coords,
                charge,
                multiplicity,
                node_counter.idx,
                SOLVENT,
                n_cores,
            )
            results.append(f"{smi_str}: Submitted job {job.job_id}")
        else:
            result = preopt(
                atoms,
                coords,
                charge,
                multiplicity,
                node_counter.idx,
                SOLVENT,
                n_cores,
            )
            results.append(f"{smi_str}: {result}")

    for r in results:
        print(r)


# # def conf_search(species):
# #     # conf search
# #     goat = orca_calculate(
# #         atoms,
# #         preopt["opt_coords"],
# #         charge,
# #         multiplicity,
# #         options={"goat": None, "XTB2": None} | ({"alpb": SOLVENT} if SOLVENT else {}),
# #         n_cores=n_cores,
# #         memory=memory,
# #     )
# #     results["goat"] = goat
# #     # SP on ensemble
# #     ensemble = goat["goat"]["ensemble"]
# #     atoms = ensemble["atoms"]
# #     sp_energies = []
# #     for coords in ensemble["coords"]:
# #         sp = orca_calculate(
# #             atoms,
# #             coords,
# #             charge,
# #             multiplicity,
# #             options=(
# #                 ENSEMBLE_SP_OPTIONS | {"SMD": solvent}
# #                 if solvent
# #                 else ENSEMBLE_SP_OPTIONS
# #             ),
# #             n_cores=n_cores,
# #             memory=memory,
# #         )
# #         if sp["normal_termination"]:
# #             sp_energies.append(sp["electronic_energy"])
# #         else:
# #             sp_energies.append(np.inf)
# #     results["ensemble_sp"] = sp_energies
# #     min_idx = np.argmin(sp_energies)
# #     coords = ensemble["coords"][min_idx]


def set_md_options():
    options = [
        MDOptions(time=10, shake=0, temp=TEMPERATURE),
        MetaDynOptions(kpush=0.15, alp=0.3),
        WallOptions(),
        SCCOptions(temp=9000),
    ]
    DEFAULT_OPTIONS = "\n".join([str(o) for o in options]) + "\n$cma"
    if typer.confirm(
        f"\nDo you want to change the MD options?\n\n{DEFAULT_OPTIONS}\n\n",
        default=False,
    ):
        PROMT = "# Options for MD:"
        message = click.edit(PROMT + "\n" + DEFAULT_OPTIONS)
        if message is not None:
            OPTIONS = "\n".join(
                [line for line in message.splitlines() if not line.startswith("#")]
            )
        OPTIONS = DEFAULT_OPTIONS
    else:
        OPTIONS = DEFAULT_OPTIONS

    return OPTIONS


def split_complex(
    atoms, coords, overall_charge, xtb_options={"alpb": "water"}, n_cores=1, scr="."
):
    def fragment_charge(total_charge, n_fragments, min_max_frag_charge=(-2, 2)):
        min_max_frag_charge = list(min_max_frag_charge)
        min_max_frag_charge[1] += 1
        min_max_frag_charge = range(*tuple(min_max_frag_charge))
        return [
            combo
            for combo in product(min_max_frag_charge, repeat=n_fragments)
            if sum(combo) == total_charge
        ]

    complex = ac2mol(
        atoms,
        coords,
        charge=overall_charge,
        use_xtb=True,
    )

    frags = list(Chem.GetMolFrags(complex, asMols=True))
    possible_charges = fragment_charge(overall_charge, len(frags))
    system_energies = []
    possible_multiplicities = []
    for pc in possible_charges:
        tmp_energies = []
        tmp_multiplicities = []
        for i, frag in enumerate(frags):
            atoms = [a.GetSymbol() for a in frag.GetAtoms()]
            coords = frag.GetConformer().GetPositions()
            charge = pc[i]
            multiplicity = (
                sum([a.GetAtomicNum() for a in frag.GetAtoms()]) - charge
            ) % 2 + 1
            tmp_multiplicities.append(multiplicity)
            results = xtb_calculate(
                atoms,
                coords,
                charge=charge,
                multiplicity=multiplicity,
                options=xtb_options,
                n_cores=n_cores,
                scr=scr,
            )
            if results["normal_termination"]:
                energy = results["electronic_energy"]
            else:
                energy = math.inf
            tmp_energies.append(energy)
        system_energies.append(sum(tmp_energies))
        possible_multiplicities.append(tmp_multiplicities)
    system_energies = np.asarray(system_energies)
    system_energies -= system_energies.min()
    system_energies *= hartree2kcalmol(1)
    best_charges = possible_charges[np.argmin(system_energies)]
    best_multiplicities = possible_multiplicities[np.argmin(system_energies)]
    return frags, best_charges, best_multiplicities


def position_species(species_list):
    """Position multiple species in 3D space to maximize separation while
    minimizing volume.

    Returns combined atoms and coordinates.
    """
    # Calculate molecular sizes and centers
    mol_sizes = []
    mol_centers = []
    all_atoms = []
    all_coords = []

    for species in species_list:
        # Get atoms and coordinates
        atoms = species.atoms
        coords = np.asarray(species.coords)

        # Calculate molecular size (maximum distance between any two atoms)
        dm = np.zeros((len(atoms), len(atoms)))
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dm[i, j] = dm[j, i] = np.linalg.norm(coords[i] - coords[j])
        mol_size = np.max(dm)
        mol_sizes.append(mol_size)

        # Calculate center of mass
        center = np.mean(coords, axis=0)
        mol_centers.append(center)

        all_atoms.extend(atoms)
        all_coords.extend(coords)

    # Position molecules in a sphere
    n_species = len(species_list)
    if n_species > 1:
        # Calculate optimal radius based on molecular sizes
        avg_size = np.mean(mol_sizes)
        radius = avg_size * (
            n_species ** (1 / 3)
        )  # Scale with cube root of number of species

        # Generate points on a sphere using the Fibonacci sphere algorithm
        points = []
        phi = np.pi * (3.0 - np.sqrt(5.0))  # Golden angle in radians
        for i in range(n_species):
            y = 1 - (i / float(n_species - 1)) * 2  # y goes from 1 to -1
            radius_i = np.sqrt(1 - y * y)  # radius at y
            theta = phi * i  # golden angle increment
            x = np.cos(theta) * radius_i
            z = np.sin(theta) * radius_i
            points.append([x, y, z])

        # Scale points to desired radius and translate molecules
        for i in range(n_species):
            if i > 0:  # Keep first molecule at origin
                translation = np.array(points[i]) * radius
                # Update coordinates for this species
                start_idx = sum(len(s.atoms) for s in species_list[:i])
                end_idx = start_idx + len(species_list[i].atoms)
                for j in range(start_idx, end_idx):
                    all_coords[j] += translation

    return all_atoms, np.array(all_coords)


@app.command()
def run_md(
    species: list[Path] = typer.Argument(None, help="List of file paths"),
    multiplicity: int = typer.Argument(None, help="Multiplicity"),
    n_runs: int = 1,
    n_cores: int = 1,
    remote: bool = False,
):
    if not species:
        species = choose_node(node_data, multiple=True, types=["species"])
    # convert node-files to species
    species = [Node.load(node_data / f) for f in species]
    charge = sum([s.data[0].charge for s in species])
    # if multiplicity not provided, prompt the user
    if multiplicity is None:
        multiplicity = int(Prompt.ask("Multiplicity", default=1))

    detailed_str = set_md_options()
    options = {} | ({"alpb": SOLVENT} if SOLVENT else {})
    if True:
        options["etemp"] = 6000

    # Position species in 3D space
    atoms, coords = position_species(species)
    ids = ",".join([str(s.idx) for s in species])

    # here is the steo where it get submitted a bunc of times, maybe make a submit flag on the method so i can run it locally or have it run on the cluster
    def run(
        atoms, coords, charge, multiplicity, options, detailed_str, n_cores, unique_id
    ):
        _logger = logging.getLogger("tooltoad.mdd")
        _logger.setLevel(logging.INFO)
        _logger.addHandler(logging.StreamHandler())
        _logger = logging.getLogger("mdd")
        _logger.setLevel(logging.INFO)
        _logger.addHandler(logging.StreamHandler())
        SCR = os.getenv("SCRATCH", ".")
        SCR = "."
        product = md_step(
            atoms,
            coords,
            charge,
            multiplicity,
            options=options,
            detailed_input_str=detailed_str,
            n_md_cores=n_cores,
            max_products=1,
            scr=SCR,
            save_traj=True,
        )
        # only keep first product, TODO: potentially keep first N products later on
        # if len(products) > 1:
        #     products = [products[0]]
        # print("MD finished")
        # for frame, product in products:
        #     product["frame"] = frame
        #     product["traj"] = traj
        #     with open(
        #         product_data / f"{ids}_product_{frame}_{unique_id}.json", "w"
        #     ) as f:
        #         json.dump(product, f)
        with open(
            product_data / f"{ids}_product_{product['frame']}_{unique_id}.json", "w"
        ) as f:
            for k, v in product.items():
                if isinstance(v, np.ndarray):
                    product[k] = v.tolist()
            json.dump(product, f)

    if remote:
        executor.update_parameters(
            name="md",
            cpus_per_task=n_cores + 1,
        )
        if n_runs > 1:
            jobs = executor.map_array(
                run,
                *[
                    [x] * n_runs
                    for x in [
                        atoms,
                        coords,
                        charge,
                        multiplicity,
                        options,
                        detailed_str,
                        n_cores,
                    ]
                ]
                + [list(range(n_runs))],
            )
            job = jobs[0]
        else:
            job = executor.submit(
                run,
                atoms,
                coords,
                charge,
                multiplicity,
                options,
                detailed_str,
                n_cores,
                0,
            )
        result = f"Submitted job {job.job_id}"
    else:
        if n_runs > 1:
            raise NotImplementedError
        result = run(
            atoms, coords, charge, multiplicity, options, detailed_str, n_cores, 0
        )
    return result


# # start with the md method
# # make node object
# # that will dictate how the qmresults pbject should look like
# # how do i deal with additional dft data?

# # do i consider stereo info?
# # no two nodes correspond to the the same molecular species.
# # differnt types of node, complex nodes, light grey, ts nodes, in teal and species nodes in blue
# # species nodes connect to complex nodes and these go to ts nodes, complex nodes again and then species nodes.
# # each node has a flag is a conf search has been run on it
# # each node has a list of qm data objects
# # qm data objects have a dict of dicts
# # key is the method, value is a dict from the qm_calculate, for now always orca
# # the list of qm results in the node is always sorted based on energy, probs r2scan-3c/gfn2xtb gibbs energy
# # control has a method where i can select all the nodes that are missing higher lvl data and submit the calcs automatically
# # do i have a way of opening the show via ssh from steno

# # TODO: am here rn

"""
reaction have reactantids -- productids
reactant and product ids can be concatenated reactant/product ind
possible separators: - _ # {}
use , to seperate id in ids, should only be permanently necessary in reactions and complexes
when its a multi, reactant or multi product reaction, a complex node is created for it and the reaction is from complex to complex or species etc.
then there must be a new function that collects all complexes and splits them. them a split is done, the reaction must be appended with the new species ids
also, when mu,tiple input species are selected, a complex node must be created. it should be created from a geometry before the reaction occurs. how do i do that, that would need to be in the run-md method
okay so i only construct the reactant complex from the ts as relaxed endpoint from irc

therefore:
1. make a complex node from reaction is it containes more then one fragment
2. in collect products, if its a complex, also create the associated species nodes and add the info the the relevant reactions
3. still, i should get the trajectory in order to extract a good alignment for the fragments


"""


# # control has a method that combs through all the new product data and adds it to the right nodes and increment the reaction counter accordigly, if no matching node exists it will make a new node
# # control has a method to generate new proton transfered nodes and also ring based H* abstraction nodes
# # control also has a method that takes a reactant complex node and makes species nodes out of it, including conf search and all
# # control has a method to run conf search on a node
# # control has a method to run a qm calculation on a node, low priority
# # control has a synch method that pulls from steno
# # each node in the nodes dir has a unique id
# # a reaction has a type, uu, ub, bu, bb, etc.
# # each reaction has reactants as a list, this is the id of the species nodes
# # also has products list, this is the id of the species nodes
# # also has a ts node, this is the id of the ts node
# # if b or t, then it also has reactant complex and product complex which is always just one id
# #
# # reaction has a IRC method to verify reaction
# # reaction has a plot method for the profile
# # reaction has a ts view_method
# # reaction has a counter on how many times it has been found

# # types of reactions:
# # consittues, from species to complex, no ts
# # hardcoded, from whatever to whatever without ts
# # reaction, has a ts, from whatever to whatever

# # the comb method could first check for n atoms, then chemical formula, then same molecular graph to make it quick
# # it only is a proper match if the graphs are isomorphic
# # node has a property which is how often it was found, this must be incremented when the reaction counter is incremented, so maybe it just pulls from the counter of all reactions ending there.
# # how do i archive that reactions are directional?
# # nodes should has a set of reaction ids that end there and that start there, maybe just pulling form the reactions dir, from the filenames
# # maybe reaction is just one json file, and its a list of dicts where each dict is a reaction with type, etc.


class Reaction:
    def __init__(
        self,
        reactant_id,
        product_id,
        reactant_species=None,
        product_species=None,
        count=1,
        rxn_type: str = "reaction",
    ):
        self.reactant_id = reactant_id
        self.product_id = product_id
        self.reactant_species = reactant_species
        self.product_species = product_species
        self.id = f"{reactant_id}-{product_id}"
        self.count = count
        self.rxn_type = rxn_type

    def get_energy_profile(self):
        # reaction has a method to get the energy of the reactants and products, this is the sum of the species nodes
        pass

    def find_ts(self):
        pass

    def run_irc(self):
        pass

    def plot(self):
        pass

    def show_ts(self):
        pass

    def __dict__(self):
        return {
            "reactant_id": self.reactant_id,
            "product_id": self.product_id,
            "count": self.count,
            "product_species": self.product_species,
            "reactant_species": self.reactant_species,
            "rxn_type": self.rxn_type,
        }

    @classmethod
    def from_dict(cls, data):
        return cls(
            reactant_id=data["reactant_id"],
            product_id=data["product_id"],
            count=data["count"],
            reactant_species=data.get("reactant_species"),
            product_species=data.get("product_species"),
            rxn_type=data.get("rxn_type", "reaction"),
        )

    def save(self, dir: Path):
        filename = f"{self.id}.json"
        with open(dir / filename, "w") as f:
            json.dump(self.__dict__(), f)
        return filename

    @classmethod
    def load(cls, filename: str):
        with open(filename, "r") as f:
            # data = json.load(f)
            decoder = json.JSONDecoder()
            content = f.read()
            content = remove_multiline_json_keys(content, "json")
            data, _ = decoder.raw_decode(content)
        return cls.from_dict(data)


def remove_multiline_json_keys(text, key="json"):
    pattern = rf'"{re.escape(key)}"\s*:\s*(\{{|\[)'  # match "json": { or "json": [
    result = []
    i = 0

    while i < len(text):
        match = re.search(pattern, text[i:])
        if not match:
            result.append(text[i:])
            break

        start = i + match.start()
        open_char = match.group(1)
        close_char = "}" if open_char == "{" else "]"

        result.append(text[i:start])  # keep everything before "json":

        # Find the matching closing brace/bracket
        j = start + match.end() - match.start()
        depth = 1
        while j < len(text) and depth > 0:
            if text[j] == open_char:
                depth += 1
            elif text[j] == close_char:
                depth -= 1
            j += 1

        # Skip trailing whitespace and comma after the block
        while j < len(text) and text[j] in " \t\r\n":
            j += 1
        if j < len(text) and text[j] == ",":
            j += 1

        i = j  # continue after removed block

    return "".join(result)


if __name__ == "__main__":
    app()
