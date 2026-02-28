import json
import os
from pathlib import Path

import numpy as np
from rdkit import Chem
from rxnnet.orca import orca_calculate

LEVELS = {
    "quick": {"opt": "XTB2", "sp": "R2SCAN-3C"},
    "normal": {"opt": "R2SCAN-3C", "sp": "wB97X-3C"},
}

OPT_SOLVENT_MODEL = {"quick": "alpb", "normal": "smd"}

N_CORES = 4
MEMORY_PER_CORE = 2  # GB

# TODO: add conf search


def optimize(
    name,
    atoms,
    coords,
    charge,
    multiplicity,
    level,
    data_dir,
    solvent=None,
    n_cores=1,
    memory=4,
):
    data_dir = Path(data_dir)
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
        # write the json file to the data directory
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
    data_dir,
    solvent=None,
    n_cores=1,
    memory=4,
):
    data_dir = Path(data_dir)
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
        # write the json file to the data directory
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
    data_dir,
    solvent=None,
    n_cores=1,
    memory=4,
):
    data_dir = Path(data_dir)
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
        # write the json file to the data directory
        json_file = data_dir / (str(name) + f"-sp-{level}.json")
        with open(json_file, "w") as f:
            json.dump(results["json"], f, indent=4)
    return results


def run_calcs(sdf_str, name, multiplicity, level, solvent, data_dir, n_cores, memory):
    data_dir = Path(data_dir)
    mol = Chem.MolFromMolBlock(sdf_str, removeHs=False)
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = mol.GetConformer().GetPositions()
    charge = Chem.GetFormalCharge(mol)

    # check if opt result already exists
    opt_file = data_dir / (str(name) + f"-opt-{level}.json")
    if opt_file.exists():
        print("")
        with open(opt_file, "r") as f:
            data = json.load(f)
        final_geom = data["Geometries"][-1]
        xyz_lines = final_geom["Geometry"]["Coordinates"]["Cartesians"]
        atoms = [line[0] for line in xyz_lines]
        opt_coords = np.array([line[1:] for line in xyz_lines]) * 0.529177
    else:
        opt_results = optimize(
            name=name,
            atoms=atoms,
            coords=coords,
            charge=charge,
            multiplicity=multiplicity,
            level=level,
            solvent=solvent,
            data_dir=data_dir,
            n_cores=n_cores,
            memory=memory,
        )
        final_geom = opt_results["json"]["Geometries"][-1]
        xyz_lines = final_geom["Geometry"]["Coordinates"]["Cartesians"]
        atoms = [line[0] for line in xyz_lines]
        opt_coords = np.array([line[1:] for line in xyz_lines]) * 0.529177

    # check if freq file is missing
    freq_file = data_dir / (str(name) + f"-freq-{level}.json")
    if not freq_file.exists():
        _ = frequencies(
            name=name,
            atoms=atoms,
            coords=opt_coords,
            charge=charge,
            multiplicity=multiplicity,
            level=level,
            solvent=solvent,
            data_dir=data_dir,
            n_cores=n_cores,
            memory=memory,
        )

    # check if sp file is missing
    freq_file = data_dir / (str(name) + f"-sp-{level}.json")
    if not freq_file.exists():
        _ = singlepoint(
            name=name,
            atoms=atoms,
            coords=opt_coords,
            charge=charge,
            multiplicity=multiplicity,
            level=level,
            solvent=solvent,
            data_dir=data_dir,
            n_cores=n_cores,
            memory=memory,
        )


if __name__ == "__main__":
    from rxnnet.config import Config
    import argparse

    parser = argparse.ArgumentParser(
        description="Perform QM calculations for nodes in network"
    )
    parser.add_argument(
        "node_files",
        nargs="*",
        help="Path(s) to node JSON file(s). If not provided, interactive selection.",
    )
    parser.add_argument(
        "-d",
        "--network-dir",
        default=".",
        help="Network directory path (default: current directory)",
    )
    parser.add_argument(
        "-n",
        "--n-cores",
        type=int,
        default=1,
        help="Number of cores to use (default: 1)",
    )
    parser.add_argument(
        "-m",
        "--memory",
        type=int,
        default=4,
        help="Memory (default: 4 GB)",
    )
    parser.add_argument(
        "-l",
        "--level",
        type=str,
        default="quick",
        help="Level of theory to use (default: quick)",
    )

    args = parser.parse_args()

    config = Config(args.network_dir)

    for node_file in args.node_files:
        # check if all three files already exist
        if all(
            [
                (
                    config.qm_data
                    / (node_file.split("-")[0] + f"-{calc_type}-{args.level}.json")
                ).exists()
                for calc_type in ["opt", "freq", "sp"]
            ]
        ):
            print(f"All calculations are finished for {node_file}")
            continue
        with open(node_file, "r") as f:
            data = json.load(f)

        sdf_str = data["data"]["gfn2-xtb"]
        name = str(data["idx"])
        multiplicity = 1
        _ = run_calcs(
            sdf_str,
            name,
            multiplicity,
            args.level,
            config.settings["qm"]["solvent"],
            config.qm_data,
            args.n_cores,
            args.memory,
        )
