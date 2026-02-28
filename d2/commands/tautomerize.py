"""Tautomerize command for generating tautomers."""

import typer

from d2.config import NetworkConfig
from d2.network import Node
from d2.tautomers import TautomerGenerator
from d2.utils import generate_and_save, keep_strict, select_and_filter_nodes

app = typer.Typer()


@app.command("tautomerize")
def tautomerize_command(
    node_ids: list[int] = typer.Argument(
        None,
        help="Node ID(s) to process. Provide one or more IDs. If no IDs are provided, all nodes will be considered.",
    ),
    n_cores: int = typer.Option(1, help="Number of cores to use for the calculation"),
    remote: bool = typer.Option(False, help="Submit QM calculation to remote executer"),
    strict: bool = typer.Option(True, help="Use strict tautomer generation"),
):
    """Generate tautomers for selected nodes."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    def strict_filtering(node_files):
        if strict:
            filtered = keep_strict(node_files, config.reaction_data)
            print(f"kept {len(filtered)} nodes based on existing reactions")
            return filtered
        return node_files

    node_files = select_and_filter_nodes(
        config=config,
        node_ids=node_ids,
        types=["species"],
        track_file=".track_tautomers.csv",
        additional_filtering=strict_filtering,
    )

    if not node_files:
        print("No nodes to process.")
        return

    tautomer_generator = TautomerGenerator()

    if remote:
        executor = config.get_executor(
            slurm_job_name="tautomers", cpus_per_task=n_cores
        )
        jobs = []
        with executor.batch():
            for file in node_files:
                print(f"Generating tautomers for: {file.stem}")
                node = Node.load(file)
                mol = node.get_mol(config.lot)
                job = executor.submit(
                    generate_and_save,
                    tautomer_generator,
                    {"mol": mol, "solvent": config.solvent, "n_cores": n_cores},
                    [node.idx],
                    "Tautomer",
                    config.new_nodes,
                )
                node.save(config.node_data)
                jobs.append(job)
        print(f"Submitted {len(jobs)} jobs for tautomers generation.")

    else:
        for file in node_files:
            print(f"Generating tautomers for: {file.stem}")
            node = Node.load(file)
            mol = node.get_mol(config.lot)
            result = generate_and_save(
                tautomer_generator,
                {"mol": mol, "solvent": config.solvent, "n_cores": n_cores},
                [node.idx],
                "Tautomer",
                config.new_nodes,
            )
            print(result)
