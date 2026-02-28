"""Run MD command for molecular dynamics simulations."""

import json
import logging
import os
from itertools import islice

import click
import numpy as np
import typer
from rdkit.Chem import rdmolops
from rich.prompt import Prompt
from tooltoad.mdd import md_step
from tooltoad.xtb import MDOptions, MetaDynOptions, SCCOptions, WallOptions

from d2.config import NetworkConfig
from d2.network import Node
from d2.utils import select_and_filter_nodes

app = typer.Typer()

_logger = logging.getLogger("tooltoad.mdd")
_logger.setLevel(logging.INFO)
_logger.addHandler(logging.StreamHandler())


def position_species(species_list):
    """Position multiple species in 3D space to maximize separation."""
    mol_sizes = []
    all_atoms = []
    all_coords = []

    for species in species_list:
        mol = species.get_mol(species.list_lots()[0])
        atoms = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = np.asarray(mol.GetConformer().GetPositions())

        dm = np.zeros((len(atoms), len(atoms)))
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dm[i, j] = dm[j, i] = np.linalg.norm(coords[i] - coords[j])
        mol_sizes.append(np.max(dm))

        all_atoms.extend(atoms)
        all_coords.extend(coords)

    n_species = len(species_list)
    if n_species > 1:
        avg_size = np.mean(mol_sizes)
        radius = avg_size * (n_species ** (1 / 3))

        # Fibonacci sphere algorithm
        points = []
        phi = np.pi * (3.0 - np.sqrt(5.0))
        for i in range(n_species):
            y = 1 - (i / float(n_species - 1)) * 2
            radius_i = np.sqrt(1 - y * y)
            theta = phi * i
            x = np.cos(theta) * radius_i
            z = np.sin(theta) * radius_i
            points.append([x, y, z])

        for i in range(n_species):
            if i > 0:
                translation = np.array(points[i]) * radius
                start_idx = sum(
                    len(
                        species_list[j]
                        .get_mol(species_list[j].list_lots()[0])
                        .GetAtoms()
                    )
                    for j in range(i)
                )
                end_idx = start_idx + len(
                    species_list[i].get_mol(species_list[i].list_lots()[0]).GetAtoms()
                )
                for j in range(start_idx, end_idx):
                    all_coords[j] += translation

    return all_atoms, np.array(all_coords)


def set_md_options():
    config = NetworkConfig()

    options = [
        MDOptions(time=10, shake=0, temp=config.temperature),
        MetaDynOptions(kpush=0.15, alp=0.3),
        WallOptions(),
        SCCOptions(temp=6000),
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


def batched(it, size: int):
    """Yield lists of up to `size` items from `it` (Py3.10 compatible)."""
    itr = iter(it)
    while True:
        chunk = list(islice(itr, size))
        if not chunk:
            return
        yield chunk


def iter_md_tasks(
    main_nodes,
    additive_node,
    n_runs: int,
    *,
    config,
    multiplicity: int,
    options: dict,
    detailed_str: str,
    n_cores: int,
):
    """Yield argument tuples for `run_simulation_for_node`."""
    for main_node in main_nodes:
        species_list = [main_node] + ([additive_node] if additive_node else [])
        charge = sum(
            rdmolops.GetFormalCharge(s.get_mol(config.lot))  # type: ignore[name-defined]
            for s in species_list
        )
        for run_i in range(n_runs):
            yield (
                main_node,
                additive_node,
                charge,
                multiplicity,
                options,
                detailed_str,
                n_cores,
                f"{main_node.idx}_{run_i}",
            )


@app.command("run-md")
def run_md_command(
    node_ids: list[int] = typer.Argument(
        None,
        help="Node ID(s) to process. Provide one or more IDs. If no IDs are provided, all nodes will be considered.",
    ),
    multiplicity: int = typer.Option(None, help="Multiplicity"),
    n_runs: int = typer.Option(100, help="Number of MD runs per node"),
    n_cores: int = typer.Option(1, help="Number of cores"),
    remote: bool = typer.Option(False, help="Submit to remote executor"),
    additive: int
    | None = typer.Option(
        None, help="Node ID to add to each simulation as additional species"
    ),
):
    """Run molecular dynamics simulations on selected species."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    additive_node = None
    if additive is not None:
        additive_files = list(config.node_data.glob(f"{additive}-*.json"))
        if additive_files:
            additive_node = Node.load(additive_files[0])
            print(f"Using additive node {additive} in all simulations")
        else:
            print(f"Warning: Additive node {additive} not found")
            return

    if additive_node:
        tracking_file = f".track_md_additive_{additive}.csv"
    else:
        tracking_file = ".track_md.csv"

    node_files = select_and_filter_nodes(
        config=config,
        node_ids=node_ids,
        types=["species"],
        track_file=tracking_file,
    )

    if not node_files:
        if additive_node:
            print(
                f"All selected nodes have already been processed with additive {additive}."
            )
        else:
            print("All selected nodes have already been processed individually.")
        return

    main_nodes = [Node.load(f) for f in node_files]

    if additive_node:
        print(
            f"Processing {len(main_nodes)} nodes, each combined with additive node {additive}"
        )
    else:
        print(f"Processing {len(main_nodes)} nodes individually")

    if multiplicity is None:
        multiplicity = int(Prompt.ask("Multiplicity", default=1))

    detailed_str = set_md_options()
    options = {} | ({"alpb": config.solvent} if config.solvent else {})
    options["etemp"] = 6000

    def run_simulation_for_node(
        main_node,
        additive_node,
        charge,
        multiplicity,
        options,
        detailed_str,
        n_cores,
        run_id,
    ):
        """Run a single MD simulation for a specific node combination."""
        scr = os.getenv("SCRATCH", ".")

        species_list = [main_node]
        if additive_node:
            species_list.append(additive_node)

        atoms, coords = position_species(species_list)
        ids = ",".join([str(s.idx) for s in species_list])

        product = md_step(
            atoms,
            coords,
            charge,
            multiplicity,
            options=options,
            detailed_input_str=detailed_str,
            n_md_cores=n_cores,
            max_products=1,
            scr=scr,
            save_traj=True,
        )

        output_file = (
            config.product_data / f"{ids}_product_{product['frame']}_{run_id}.json"
        )
        with open(output_file, "w") as f:
            for k, v in product.items():
                if isinstance(v, np.ndarray):
                    product[k] = v.tolist()
            json.dump(product, f)

        return f"MD simulation completed. Product saved as {output_file.name}"

    print(f"Starting {n_runs} MD run(s) for each of {len(main_nodes)} node(s)")

    if remote:
        executor = config.get_executor(
            slurm_job_name="md",
            cpus_per_task=n_cores + 1,
        )
        max_array = 1000
        jobs = []
        total_count = 0
        task_iter = iter_md_tasks(
            main_nodes,
            additive_node,
            n_runs,
            config=config,
            multiplicity=multiplicity,
            options=options,
            detailed_str=detailed_str,
            n_cores=n_cores,
        )

        for chunk in batched(task_iter, max_array):
            with executor.batch():
                for args in chunk:
                    job = executor.submit(run_simulation_for_node, *args)  # type: ignore[name-defined]
                    jobs.append(job)
            total_count += len(chunk)

        print(
            f"Submitted {total_count} MD jobs "
            f"({n_runs} per node) as {max(1, (total_count + max_array - 1)//max_array)} arrays "
            f"(≤{max_array} tasks/array" + ")."
        )
    else:
        for main_node in main_nodes:
            species_list = [main_node]
            if additive_node:
                species_list.append(additive_node)

            charge = sum(
                [rdmolops.GetFormalCharge(s.get_mol(config.lot)) for s in species_list]
            )

            if n_runs > 1:
                print(
                    f"Multiple local runs not yet implemented. Running single simulation for node {main_node.idx}."
                )

            result = run_simulation_for_node(
                main_node,
                additive_node,
                charge,
                multiplicity,
                options,
                detailed_str,
                n_cores,
                f"{main_node.idx}_0",
            )
            print(result)
