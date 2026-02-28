"""Deprotonate command for generating deprotonated structures."""

import typer
from rdkit.Chem import rdmolops

from d2.config import NetworkConfig
from d2.network import Node
from d2.protonation import Protonator
from d2.utils import generate_and_save, keep_strict, select_and_filter_nodes

app = typer.Typer()


@app.command("deprotonate")
def deprotonate_command(
    node_ids: list[int] = typer.Argument(
        None,
        help="Node ID(s) to process. Provide one or more IDs. If no IDs are provided, all nodes will be considered.",
    ),
    n_cores: int = typer.Option(1, help="Number of cores to use for the calculation"),
    remote: bool = typer.Option(False, help="Submit QM calculation to remote executer"),
    min_charge: int = typer.Option(-1, help="Minimum charge to consider"),
    strict: bool = typer.Option(True, help="Use strict deprotonation generation"),
):
    """Generate deprotonated structures for selected nodes."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    def strict_filtering(node_files):
        if strict:
            return keep_strict(node_files, config.reaction_data)
        return node_files

    node_files = select_and_filter_nodes(
        config=config,
        node_ids=node_ids,
        types=["species"],
        track_file=".track_deprotonation.csv",
        additional_filtering=strict_filtering,
    )

    if not node_files:
        print("No nodes to process.")
        return

    deprotonator = Protonator(mode="deprotonate")

    if remote:
        executor = config.get_executor(
            slurm_job_name="deprotonation", cpus_per_task=n_cores
        )
        jobs = []
        with executor.batch():
            for file in node_files:
                print(f"Generating Deprotonated Structure for: {file.stem}")
                node = Node.load(file)
                mol = node.get_mol(config.lot)
                charge = rdmolops.GetFormalCharge(mol)
                if charge <= min_charge:
                    print(
                        f"Skipping {file.stem} with charge {charge} (min allowed: {min_charge})"
                    )
                    continue
                job = executor.submit(
                    generate_and_save,
                    deprotonator,
                    {"mol": mol, "solvent": config.solvent, "n_cores": n_cores},
                    [node.idx],
                    "Deprotonation",
                    config.new_nodes,
                )
                jobs.append(job)
        print(f"Submitted {len(jobs)} jobs for Deprotonation generation.")

    else:
        for file in node_files:
            print(f"Generating Deprotonated Structure for: {file.stem}")
            node = Node.load(file)
            mol = node.get_mol(config.lot)
            charge = rdmolops.GetFormalCharge(mol)
            if charge <= min_charge:
                print(
                    f"Skipping {file.stem} with charge {charge} (min allowed: {min_charge})"
                )
                continue
            result = generate_and_save(
                deprotonator,
                {"mol": mol, "solvent": config.solvent, "n_cores": n_cores},
                [node.idx],
                "Deprotonation",
                config.new_nodes,
            )
            print(result)
