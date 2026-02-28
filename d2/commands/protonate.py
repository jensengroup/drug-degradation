"""Protonate command for generating protonated structures."""

import typer
from rdkit.Chem import rdmolops

from d2.config import NetworkConfig
from d2.network import Node
from d2.protonation import Protonator
from d2.utils import generate_and_save, keep_strict, select_and_filter_nodes

app = typer.Typer()


@app.command("protonate")
def protonate_command(
    node_ids: list[int] = typer.Argument(
        None,
        help="Node ID(s) to process. Provide one or more IDs. If no IDs are provided, all nodes will be considered.",
    ),
    n_cores: int = typer.Option(1, help="Number of cores to use for the calculation"),
    remote: bool = typer.Option(False, help="Submit QM calculation to remote executer"),
    max_charge: int = typer.Option(1, help="Maximum charge to consider"),
    strict: bool = typer.Option(True, help="Use strict protonation generation"),
):
    """Generate protonated structures for selected nodes."""
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
        track_file=".track_protonation.csv",
        additional_filtering=strict_filtering,
    )

    if not node_files:
        print("No nodes to process.")
        return

    protonator = Protonator(mode="protonate")

    if remote:
        executor = config.get_executor(
            slurm_job_name="protonation", cpus_per_task=n_cores
        )
        jobs = []
        with executor.batch():
            for file in node_files:
                print(f"Generating Protonated Structure for: {file.stem}")
                node = Node.load(file)
                mol = node.get_mol(config.lot)
                charge = rdmolops.GetFormalCharge(mol)
                if charge >= max_charge:
                    print(
                        f"Skipping {file.stem} with charge {charge} (max allowed: {max_charge})"
                    )
                    continue
                job = executor.submit(
                    generate_and_save,
                    protonator,
                    {"mol": mol, "solvent": config.solvent, "n_cores": n_cores},
                    [node.idx],
                    "Protonation",
                    config.new_nodes,
                )
                jobs.append(job)
        print(f"Submitted {len(jobs)} jobs for Protonation generation.")

    else:
        for file in node_files:
            print(f"Generating Protonated Structure for: {file.stem}")
            node = Node.load(file)
            mol = node.get_mol(config.lot)
            charge = rdmolops.GetFormalCharge(mol)
            if charge >= max_charge:
                print(
                    f"Skipping {file.stem} with charge {charge} (max allowed: {max_charge})"
                )
                continue
            result = generate_and_save(
                protonator,
                {"mol": mol, "solvent": config.solvent, "n_cores": n_cores},
                [node.idx],
                "Protonation",
                config.new_nodes,
            )
            print(result)
