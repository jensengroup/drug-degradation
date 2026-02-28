"""Add command for adding new molecules to the network."""

import os
from typing import List

import typer
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdmolops
from tooltoad.orca import orca_calculate

from d2.config import IdxCounter, NetworkConfig
from d2.network import Node
from d2.utils import make_mol

app = typer.Typer()


def preopt(atoms, coords, charge, multiplicity, idx, solvent, n_cores, config):
    """Pre-optimization function for new molecules."""
    scr = os.getenv("SCRATCH", ".")
    results = orca_calculate(
        atoms,
        coords,
        charge,
        multiplicity,
        options={"opt": None, "XTBFF": None} | ({"alpb": solvent} if solvent else {}),
        scr=scr,
        n_cores=n_cores,
    )
    results = orca_calculate(
        atoms,
        results["opt_coords"],
        charge,
        multiplicity,
        options={"opt": None, "XTB2": None} | ({"alpb": solvent} if solvent else {}),
        scr=scr,
        n_cores=n_cores,
    )
    # make species
    mol = make_mol(results, "opt_coords")
    mol.SetProp("origin-ids", str([0]))
    mol.SetProp("origin-type", "add")
    node = Node(
        data={"gfn2-xtb": mol},
        idx=idx,
    )
    node.save(config.node_data)
    return mol


@app.command("add")
def add_command(
    smi: List[str] = typer.Argument(
        ..., help="One or more SMILES strings for molecules"
    ),
    multiplicity: int = typer.Option(1, help="Multiplicity of the Molecule"),
    remote: bool = typer.Option(False, help="Submit QM calculation to remote executer"),
    n_cores: int = typer.Option(1, help="Number of cores to use"),
):
    """Add new molecules to the network from SMILES strings."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    node_counter = IdxCounter(str(config.network_file), "node_count")

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
            executor = config.get_executor(
                slurm_job_name="preopt-species", cpus_per_task=n_cores
            )
            job = executor.submit(
                preopt,
                atoms,
                coords,
                charge,
                multiplicity,
                node_counter.idx,
                config.solvent,
                n_cores,
                config,
            )
            results.append(f"{smi_str}: Submitted job {job.job_id}")
        else:
            try:
                mol = preopt(
                    atoms,
                    coords,
                    charge,
                    multiplicity,
                    node_counter.idx,
                    config.solvent,
                    n_cores,
                    config,
                )
                results.append(
                    f"{smi_str}: Successfully added as node {node_counter.idx}"
                )
            except Exception as e:
                results.append(f"{smi_str}: Failed - {e}")

    for r in results:
        print(r)
