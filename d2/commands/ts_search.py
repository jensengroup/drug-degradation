"""Transition state search command."""

import json
import os
from pathlib import Path

import click
import numpy as np
import typer
from rdkit import Chem
from rdkit.Chem import rdmolops
from tooltoad.chemutils import (
    Constraint,
    ac2mol,
    get_bond_change,
    get_connectivity_smiles,
)
from tooltoad.ndscan import PotentialEnergySurface, ScanCoord, orca_calculate
from tooltoad.ts_utils import (
    get_scan_ts_guess,
    get_ssm_ts_guess,
    sort_start_end,
)

from d2.config import NetworkConfig


def irc_check(forward, backward, reactant, product):
    return set([get_connectivity_smiles(mol) for mol in [forward, backward]]) == set(
        [get_connectivity_smiles(mol) for mol in [reactant, product]]
    )


app = typer.Typer()


def is_relevant(reactant_file, product_file):
    with open(reactant_file, "r") as f:
        reactant_data = json.load(f)
    with open(product_file, "r") as f:
        product_data = json.load(f)
    reactant = Chem.MolFromMolBlock(reactant_data["data"]["gfn2-xtb"], removeHs=False)
    product = Chem.MolFromMolBlock(product_data["data"]["gfn2-xtb"], removeHs=False)
    if reactant.GetNumAtoms() != product.GetNumAtoms():
        return False
    bond_changes = get_bond_change(reactant, product)
    atoms = [a.GetSymbol() for a in product.GetAtoms()]
    if all([any([atoms[idx] == "H" for idx in bc[1]]) for bc in bond_changes]):
        return False
    return True


def ts_search(
    reactant_file, product_file, n_cores=4, qm_data_dir=".", method="auto", reverse=True
):
    scr = os.getenv("SCRATCH", ".")
    with open(reactant_file, "r") as f:
        reactant_data = json.load(f)
    with open(product_file, "r") as f:
        product_data = json.load(f)
    reactant = Chem.MolFromMolBlock(reactant_data["data"]["gfn2-xtb"], removeHs=False)
    product = Chem.MolFromMolBlock(product_data["data"]["gfn2-xtb"], removeHs=False)
    charge = rdmolops.GetFormalCharge(product)
    start, end = sort_start_end(reactant, product, reverse=reverse)
    bond_changes = get_bond_change(start, end)
    print(f"Reaction between {start} and {end}")
    print(f"Bond Changes: {bond_changes}")
    atoms = [a.GetSymbol() for a in start.GetAtoms()]
    coords = start.GetConformer().GetPositions()
    success = False
    assert (
        len(bond_changes) > 0
    ), "No bond changes detected between reactant and product."
    if len(bond_changes) == 1 and method.lower() != "ssm":
        print("Creating TS guess with coordinate scan")
        ts_guess = get_scan_ts_guess(
            atoms,
            coords,
            bond_changes,
            charge=charge,
            xtb_options={"alpb": "water"},
            scr=scr,
        )
    else:
        print("Creating TS guess with SSM")
        ts_guess = get_ssm_ts_guess(
            atoms,
            coords,
            bond_changes,
            charge=charge,
            orca_options={"XTB2": None, "alpb": "water"},
            gsm_executable="~/opt/GSM/gsm.orca",
            scr=scr,
        )

    input_str = """%geom
    Calc_Hess true
    end"""
    ts_atoms = [a.GetSymbol() for a in ts_guess.GetAtoms()]
    ts_coords = ts_guess.GetConformer().GetPositions()

    print("TS optimization...")
    ts_data = orca_calculate(
        ts_atoms,
        ts_coords,
        charge=charge,
        options={"optts": None, "freq": None, "xtb2": None, "alpb": "water"},
        xtra_inp_str=input_str,
        n_cores=n_cores,
        memory=4 * n_cores,
        scr=scr,
    )

    print("Run IRC...")
    irc = orca_calculate(
        ts_atoms,
        ts_data["opt_coords"],
        charge=charge,
        options={"xtb2": None, "alpb": "water", "irc": None},
        n_cores=n_cores,
        memory=4 * n_cores,
        scr=scr,
    )
    forward = ac2mol(
        irc["irc"]["forward"]["atoms"], irc["irc"]["forward"]["opt_coords"]
    )
    backward = ac2mol(
        irc["irc"]["backward"]["atoms"], irc["irc"]["backward"]["opt_coords"]
    )

    if irc_check(forward, backward, reactant, product):
        success = True
    else:
        print("IRC did not connect reactant and product")
    if not success:
        new_endpoints = []
        for endpoint in [forward, backward]:
            if get_connectivity_smiles(endpoint) not in [
                get_connectivity_smiles(mol) for mol in [reactant, product]
            ]:
                print("Optimizing Endpoint...")
                atoms = [a.GetSymbol() for a in endpoint.GetAtoms()]
                coords = endpoint.GetConformer().GetPositions()
                results = orca_calculate(
                    atoms,
                    coords,
                    charge=charge,
                    options={"opt": None, "xtb2": None, "alpb": "water"},
                    n_cores=n_cores,
                    memory=4 * n_cores,
                    scr=scr,
                )
                endpoint = ac2mol(atoms, results["opt_coords"])
                new_endpoints.append(endpoint)
            else:
                new_endpoints.append(endpoint)
        if irc_check(*new_endpoints, reactant, product):
            success = True
        else:
            with open("fail_log.txt", "a") as f:
                f.write(f"{start},{end}:{bond_changes}\n")
            raise RuntimeError("Multi Step reaction, need manual intervention")
    if not success:
        raise RuntimeError("All attempts to find the TS have failed")
    print("Found the right TS")
    with open(
        Path(qm_data_dir)
        / f"ts-{reactant_file.stem.split('-')[0]}-{product_file.stem.split('-')[0]}.json",
        "w",
    ) as f:
        for k, v in ts_data.items():
            if isinstance(v, np.ndarray):
                ts_data[k] = v.tolist()
        json.dump(ts_data, f)
    # to r2scan-3c sp on ts
    sp = orca_calculate(
        atoms,
        ts_data["opt_coords"],
        charge=charge,
        options={"r2scan-3C": None, "smd": "water"},
        n_cores=n_cores,
        memory=4 * n_cores,
    )
    # now write the ts file
    mol = ac2mol(atoms, ts_data["opt_coords"], charge=charge)
    mol.SetDoubleProp("l1_gibbs-energy", ts_data["gibbs_energy"])
    mol.SetDoubleProp("l1_electronic-energy", ts_data["electronic_energy"])
    mol.SetDoubleProp("l2_electronic-energy", sp["electronic_energy"])
    mol.SetDoubleProp(
        "l2l1_gibbs-energy",
        sp["electronic_energy"]
        + ts_data["gibbs_energy"]
        - ts_data["electronic_energy"],
    )
    qm_data_dir = Path(qm_data_dir)
    with Chem.SDWriter(
        qm_data_dir
        / f"ts-{reactant_file.stem.split('-')[0]}-{product_file.stem.split('-')[0]}.sdf"
    ) as writer:
        writer.write(mol)

    return mol


def ts_search_v1(reactant_file, product_file, n_cores=4, qm_data_dir="."):
    with open(reactant_file, "r") as f:
        reactant_data = json.load(f)
    with open(product_file, "r") as f:
        product_data = json.load(f)
    reactant = Chem.MolFromMolBlock(reactant_data["data"]["gfn2-xtb"], removeHs=False)
    product = Chem.MolFromMolBlock(product_data["data"]["gfn2-xtb"], removeHs=False)
    charge = rdmolops.GetFormalCharge(product)
    # cleanup stereo and bond order for irc check
    reactant = ac2mol(
        [a.GetSymbol() for a in reactant.GetAtoms()],
        reactant.GetConformer().GetPositions(),
    )
    product = ac2mol(
        [a.GetSymbol() for a in product.GetAtoms()],
        product.GetConformer().GetPositions(),
    )
    # get bond diff
    ac1 = rdmolops.GetAdjacencyMatrix(reactant)
    ac2 = rdmolops.GetAdjacencyMatrix(product)
    diff = ac2 - ac1

    if np.abs(diff).sum() == 2:
        pass
    else:
        raise NotImplementedError

    aIds = np.where(diff != 0)[0]
    bond_change = diff[aIds[0], aIds[1]]

    atoms = [a.GetSymbol() for a in product.GetAtoms()]
    coords = product.GetConformer().GetPositions()

    preopt = orca_calculate(
        atoms,
        coords,
        charge,
        options={"XTB2": None, "alpb": "water", "tightopt": None},
        n_cores=n_cores,
        memory=4 * n_cores,
    )
    opt_coords = preopt["opt_coords"]

    scs = [
        ScanCoord.from_current_position(
            atoms,
            opt_coords,
            aIds,
            50,
            bool(bond_change),
        )
    ]

    pes = PotentialEnergySurface(atoms, opt_coords, charge, scan_coords=scs)
    pes.xtb(n_cores=n_cores, xtb_options={"alpb": "water"}, max_cycle=25)
    pes.refine(n_cores=n_cores, orca_options={"r2scan-3c": None, "smd": "water"})
    ts_guess = pes.traj_tensor[pes.refined_pes_tensor.argmax()]

    input_str = f"""%geom
        Constraints
        {Constraint([int(i) for i in aIds], None).orca} end
    end
    """
    ts_preopt = orca_calculate(
        atoms,
        ts_guess,
        charge=charge,
        options={
            "r2scan-3c": None,
            "smd": "water",
            "opt": None,
        },
        xtra_inp_str=input_str,
        n_cores=n_cores,
        memory=4 * n_cores,
    )

    input_str = """%geom
    Calc_Hess true
    end"""
    ts = orca_calculate(
        atoms,
        ts_preopt["opt_coords"],
        charge=charge,
        options={"optts": None, "freq": None, "r2scan-3c": None, "smd": "water"},
        xtra_inp_str=input_str,
        n_cores=n_cores,
        memory=4 * n_cores,
    )
    #  check theres only one proper imaginary freq
    freq_check = sum(np.array([v["frequency"] for v in ts["vibs"]]) < -10) == 1
    if not freq_check:
        print([v["frequency"] for v in ts["vibs"]])
        raise ValueError
    irc = orca_calculate(
        atoms,
        ts["opt_coords"],
        charge=charge,
        options={"r2scan-3c": None, "smd": "water", "irc": None},
        n_cores=n_cores,
        memory=4 * n_cores,
    )
    forward = ac2mol(
        irc["irc"]["forward"]["atoms"], irc["irc"]["forward"]["opt_coords"]
    )
    backward = ac2mol(
        irc["irc"]["backward"]["atoms"], irc["irc"]["backward"]["opt_coords"]
    )

    irc_check = set([Chem.MolToSmiles(m) for m in [forward, backward]]) == set(
        [Chem.MolToSmiles(m) for m in [reactant, product]]
    )
    if not irc_check:
        print([Chem.MolToSmiles(m) for m in [forward, backward]])
        raise ValueError
    # to r2scan-3c sp on ts
    sp = orca_calculate(
        atoms,
        ts["opt_coords"],
        charge=charge,
        options={"wB97X-3C": None, "smd": "water"},
        n_cores=n_cores,
        memory=4 * n_cores,
    )
    # now write the ts file
    mol = ac2mol(atoms, ts["opt_coords"], charge=charge)
    mol.SetDoubleProp("l1_gibbs-energy", ts["gibbs_energy"])
    mol.SetDoubleProp("l1_electronic-energy", ts["electronic_energy"])
    mol.SetDoubleProp("l2_electronic-energy", sp["electronic_energy"])
    mol.SetDoubleProp(
        "l2l1_gibbs-energy",
        sp["electronic_energy"] + ts["gibbs_energy"] - ts["electronic_energy"],
    )
    qm_data_dir = Path(qm_data_dir)
    with Chem.SDWriter(
        qm_data_dir
        / f"ts-{reactant_file.stem.split('-')[0]}-{product_file.stem.split('-')[0]}.sdf"
    ) as writer:
        writer.write(mol)

    return mol


@app.command("ts-search")
def ts_command(
    node_ids: list[str] = typer.Argument(
        None,  # Optional argument
        help="Node ID pairs to process. Provide as '2-36 5-6 10-10000' where each pair represents reactant-product node IDs.",
    ),
    n_cores: int = typer.Option(1, help="Number of cores to use for the calculation"),
    remote: bool = typer.Option(False, help="Submit QM calculation to remote executer"),
    only_relevant: bool = typer.Option(True, help="Only consider relevant node pairs"),
    method: str = typer.Option(
        "auto", help="Method to use for TS search, either auto or ssm"
    ),
    reverse: bool = typer.Option(True, help="Whether to reverse the reaction"),
):
    """Generate transition state search for selected nodes."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    # Parse node_ids string into list of tuples
    if not node_ids:
        print("No node pairs provided. Please specify node pairs like '2-36 5-6'.")
        return

    node_id_list = []
    for pair in node_ids:
        try:
            # Split by hyphen and convert to integers
            parts = pair.split("-")
            if len(parts) != 2:
                print(
                    f"Invalid node pair format: {pair}. Expected format: 'reactant-product'"
                )
                continue
            reactant_id, product_id = int(parts[0]), int(parts[1])
            node_id_list.append((reactant_id, product_id))
        except ValueError:
            print(f"Error parsing node pair: {pair}. Both IDs must be integers.")
            continue

    # Find node files corresponding to the indices
    node_files = []
    for node_id_tuple in node_id_list:
        # Convert tuple to string format that matches file naming
        try:
            node_file_pair = [
                next(config.node_data.glob(f"{idx}-*.json")) for idx in node_id_tuple
            ]
        except StopIteration:
            print(
                f"Warning: Node file for ID {node_id_tuple} not found at {config.node_data}"
            )
            continue
        node_files.append(node_file_pair)

    if only_relevant:
        tmp = []
        for file_pair in node_files:
            if is_relevant(*file_pair):
                tmp.append(file_pair)
        node_files = tmp

    tmp = []
    for f in node_files:
        if not Path(
            config.qm_data
            / f"ts-{f[0].stem.split('-')[0]}-{f[1].stem.split('-')[0]}.sdf"
        ).exists():
            tmp.append(f)
    node_files = tmp

    if not node_files:
        print("No nodes to process.")
        return

    if remote:
        executor = config.get_executor(
            slurm_job_name="ts_search", cpus_per_task=n_cores
        )
        jobs = []
        with executor.batch():
            for pair in node_files:
                click.echo(f"TS Search for: {pair[0].stem} and {pair[1].stem}")
                job = executor.submit(
                    ts_search,
                    pair[0],
                    pair[1],
                    n_cores,
                    config.qm_data,
                    method,
                    reverse,
                )
                jobs.append(job)
        print(f"Submitted {len(jobs)} jobs for TS search generation.")

    else:
        for pair in node_files:
            click.echo(f"Generating TS Search for: {pair[0].stem} and {pair[1].stem}")
            result = ts_search(pair[0], pair[1], n_cores, config.qm_data)
            print(result)
