import json
import os
import sys
from pathlib import Path

import click
import numpy as np
import submitit
from rdkit import Chem
from tooltoad.orca import orca_calculate

LEVELS = {
    "quick": {"opt": "XTB2", "sp": "R2SCAN-3C"},
    "normal": {"opt": "R2SCAN-3C", "sp": "wB97X-3C"},
}

OPT_SOLVENT_MODEL = {"quick": "alpb", "normal": "smd"}
# SETUP

# data_dir = Path(os.getenv("DATA_DIR", "/groups/kemi/julius/data"))
# data_dir = Path(os.getenv("DATA_DIR", "/Users/julius/opt/tool-toad/scipts/data"))
# data_dir.mkdir(exist_ok=True)

N_CORES = 4
MEMORY_PER_CORE = 2  # GB

executor = submitit.AutoExecutor(folder="dft")
executor.update_parameters(
    timeout_min=6000, slurm_partition="kemi1", slurm_array_parallelism=250
)


def optimize(
    name,
    atoms,
    coords,
    charge,
    multiplicity,
    level,
    solvent=None,
    n_cores=1,
    memory=4,
    data_dir=".",
):
    scratch = os.getenv("SCRATCH", ".")
    options = {"opt": None, LEVELS[level]["opt"]: None}
    if solvent:
        options[OPT_SOLVENT_MODEL[level]] = solvent
    results = orca_calculate(
        atoms=atoms,
        coords=coords,
        charge=charge,
        multiplicity=multiplicity,
        options=options,
        scr=scratch,
        n_cores=n_cores,
        memory=memory,
        log_file=(data_dir / (str(name) + f"-opt-{level}.log")).absolute(),
    )
    if results["normal_termination"]:
        json_file = data_dir / (str(name) + f"-opt-{level}.json")
        with open(json_file, "w") as f:
            json.dump(results["json"], f, indent=4)
    return results


def frequencies(
    name,
    atoms,
    coords,
    charge,
    multiplicity,
    level,
    solvent=None,
    n_cores=1,
    memory=4,
    data_dir=".",
):
    scratch = os.getenv("SCRATCH", ".")
    options = {"freq": None, LEVELS[level]["opt"]: None}
    if solvent:
        options[OPT_SOLVENT_MODEL[level]] = solvent
    results = orca_calculate(
        atoms=atoms,
        coords=coords,
        charge=charge,
        multiplicity=multiplicity,
        options=options,
        scr=scratch,
        n_cores=n_cores,
        memory=memory,
        log_file=(data_dir / (str(name) + f"-freq-{level}.log")).absolute(),
    )
    if results["normal_termination"]:
        json_file = data_dir / (str(name) + f"-freq-{level}.json")
        with open(json_file, "w") as f:
            json.dump(results["json"], f, indent=4)
    return results


def singlepoint(
    name,
    atoms,
    coords,
    charge,
    multiplicity,
    level,
    solvent=None,
    n_cores=1,
    memory=4,
    data_dir=".",
):
    scratch = os.getenv("SCRATCH", ".")
    options = {LEVELS[level]["sp"]: None}
    if solvent:
        options["smd"] = solvent
    results = orca_calculate(
        atoms=atoms,
        coords=coords,
        charge=charge,
        multiplicity=multiplicity,
        options=options,
        scr=scratch,
        n_cores=n_cores,
        memory=memory,
        log_file=(data_dir / (str(name) + f"-sp-{level}.log")).absolute(),
    )
    if results["normal_termination"]:
        json_file = data_dir / (str(name) + f"-sp-{level}.json")
        with open(json_file, "w") as f:
            json.dump(results["json"], f, indent=4)
    return results


def optfreq(
    sdf_str, name, level, multiplicity, solvent, resource_multiplier, data_dir="."
):
    """Perform optimization, frequency, and single-point calculations."""
    data_dir = Path(data_dir)
    mol = Chem.MolFromMolBlock(sdf_str, removeHs=False)
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = mol.GetConformer().GetPositions()
    charge = Chem.GetFormalCharge(mol)
    click.echo(f"Name: {name}")
    click.echo(f"Solvent: {solvent}")

    n_cores = int(N_CORES * resource_multiplier)
    memory = int(N_CORES * MEMORY_PER_CORE * resource_multiplier)
    memory_freq = memory * 2
    if level.lower() == "normal":
        memory_freq *= 2

    opt_file = data_dir / (str(name) + f"-opt-{level}.json")
    if opt_file.exists():
        click.echo(f"Optimization result already exists: {opt_file}")
        with open(opt_file, "r") as f:
            data = json.load(f)
        final_geom = data["Geometries"][-1]
        xyz_lines = final_geom["Geometry"]["Coordinates"]["Cartesians"]
        atoms = [line[0] for line in xyz_lines]
        opt_coords = [line[1:] for line in xyz_lines]
        try:
            l1_electronic_energy = final_geom["DFT_Energy"]["finalEn"]
        except KeyError:
            l1_electronic_energy = final_geom["Single_Point_Data"]["FinalEnergy"]
    else:
        click.echo(
            f"Running geometry optimization using {n_cores} cores and {memory} GB memory."
        )
        executor.update_parameters(
            slurm_cpus_per_task=n_cores,
            slurm_mem_gb=memory,
            slurm_job_name=f"{name}-opt",
        )
        job = executor.submit(
            optimize,
            name,
            atoms,
            coords,
            charge,
            multiplicity,
            level,
            solvent,
            n_cores,
            memory,
            data_dir,
        )
        sys.exit()
        opt = job.result()
        opt_coords = opt["opt_coords"]
        l1_electronic_energy = opt["electronic_energy"]
        if not opt["normal_termination"]:
            click.echo(
                f"Optimization did not terminate normally, see log file {(data_dir / (str(name) + f'-opt-{level}.log')).absolute()}. Exiting."
            )
            return
    mol.SetProp("l1", LEVELS[level]["opt"])
    mol.GetConformer().SetPositions(np.array(opt_coords))
    mol.SetDoubleProp("l1_electronic-energy", l1_electronic_energy)

    freq_file = data_dir / (str(name) + f"-freq-{level}.json")
    if freq_file.exists():
        click.echo(f"Frequency result already exists: {freq_file}")
        with open(freq_file, "r") as f:
            data = json.load(f)
        thermo = data["Geometries"][0]["THERMOCHEMISTRY_Energies"][0]
        freq = {
            "gibbs_energy": thermo["freeEnergyG"],
            "electronic_energy": thermo["elEnergy"],
        }

    else:
        click.echo(
            f"Running frequency calculation using {n_cores} cores and {memory_freq} GB memory."
        )
        executor.update_parameters(
            slurm_cpus_per_task=n_cores,
            slurm_mem_gb=memory_freq,
            slurm_job_name=f"{name}-freq",
        )
        job = executor.submit(
            frequencies,
            name,
            atoms,
            coords,
            charge,
            multiplicity,
            level,
            solvent,
            n_cores,
            memory,
            data_dir,
        )
        sys.exit()
        freq = job.result()
        if not freq["normal_termination"]:
            click.echo(
                f"Frequency calculation did not terminate normally, see log file {(data_dir / (str(name) + f'-freq-{level}.log')).absolute()}. Exiting."
            )
            return
    mol.SetDoubleProp("l1_gibbs-energy", freq["gibbs_energy"])
    mol.SetDoubleProp(
        "l1_gibbs-correction", freq["gibbs_energy"] - freq["electronic_energy"]
    )
    executor.update_parameters(
        slurm_cpus_per_task=n_cores, slurm_mem_gb=memory, slurm_job_name=f"{name}-sp"
    )
    click.echo(
        f"Running single point calculation using {n_cores} cores and {memory} GB memory."
    )
    job = executor.submit(
        singlepoint,
        name,
        atoms,
        coords,
        charge,
        multiplicity,
        level,
        solvent,
        n_cores,
        memory,
        data_dir,
    )
    sp = job.result()
    if not sp["normal_termination"]:
        click.echo(
            f"Single point calculation did not terminate normally, see log file {(data_dir / (str(name) + f'-sp-{level}.log')).absolute()}. Exiting."
        )
        return

    mol.SetProp("l2", LEVELS[level]["sp"])
    mol.SetDoubleProp("l2_electronic-energy", sp["electronic_energy"])
    mol.SetDoubleProp(
        "l2l1_gibbs-energy",
        sp["electronic_energy"] + freq["gibbs_energy"] - freq["electronic_energy"],
    )

    with Chem.SDWriter((data_dir / (str(name) + f"-{level}.sdf")).absolute()) as writer:
        writer.write(mol)

    return mol


if __name__ == "__main__":
    optfreq()
