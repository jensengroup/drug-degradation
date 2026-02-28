import click
import numpy as np
import submitit
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdmolops
from tooltoad.chemutils import (
    ConformerCalculator,
    canonicalize_solvent,
    energy_filter_conformer,
    filter_conformers,
)
from tooltoad.orca import orca_calculate
from tooltoad.xtb import xtb_calculate

executor = submitit.AutoExecutor(folder="test")
executor.update_parameters(
    timeout_min=600, slurm_partition="kemi1", slurm_array_parallelism=10
)

Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)


def run_goat(data):
    smi = data["smi"]
    # embed initial conformer
    mol3d = Chem.AddHs(Chem.MolFromSmiles(smi))
    rdDistGeom.EmbedMolecule(mol3d, randomSeed=42)
    atoms = [a.GetSymbol() for a in mol3d.GetAtoms()]
    coords = mol3d.GetConformer().GetPositions()
    charge = rdmolops.GetFormalCharge(mol3d)
    orca_options = data.pop("orca_options")

    # preopt
    preopt_options = orca_options.copy()
    preopt_options.pop("GOAT")
    preopt_options["opt"] = None
    preopt = orca_calculate(
        atoms=atoms,
        coords=coords,
        charge=charge,
        options=preopt_options,
        n_cores=data.get("n_cores", 1),
    )

    # first get the xtb gfn2 gibbs correction
    goat = orca_calculate(
        atoms,
        preopt["opt_coords"],
        charge,
        options=orca_options,
        xtra_inp_str=data.get("detailed_input_str", ""),
        n_cores=data.get("n_cores", 1),
        read_files=["orca.out"],
    )

    return goat


def run_orca_sp(data):
    results = orca_calculate(
        data["atoms"],
        data["coords"],
        data["charge"],
        options=data["orca_options"],
        n_cores=data.get("n_cores", 1),
    )
    if results["normal_termination"]:
        energy = results["electronic_energy"]
    else:
        energy = np.inf
    return energy


def generate_with_goat(smi, solvent, n_cores):
    orca_options = {"GOAT": None, "XTB2": None}
    charge = rdmolops.GetFormalCharge(Chem.MolFromSmiles(smi))
    if solvent:
        orca_options["smd"] = canonicalize_solvent(solvent, "xtb")
    detailed_input = """
    %GOAT
    GFNUPHILL    GFNFF
    END
    """
    data = {
        "smi": smi,
        "orca_options": orca_options,
        "n_cores": n_cores,
        "detailed_input_str": detailed_input,
    }

    job = executor.submit(run_goat, data)
    click.echo("Calculations submitted. Waiting for results...")
    result = job.result()
    click.echo("Calculations completed.")

    ensemble = result["goat"]["ensemble"]
    atoms = ensemble["atoms"]
    jobs = []
    orca_sp_options = {"r2scan-3c": None}
    if solvent:
        orca_sp_options["smd"] = orca_options["smd"]
    with executor.batch():
        for c in ensemble["coords"]:
            data = {
                "atoms": atoms,
                "coords": c,
                "charge": charge,
                "orca_options": orca_sp_options,
            }
            job = executor.submit(run_orca_sp, data)
            jobs.append(job)
    click.echo("Calculations submitted. Waiting for results...")
    results = [j.result() for j in jobs]
    click.echo("Calculations completed.")
    # add sp energies to the goat dict
    result["goat"]["ensemble"]["sps"] = results
    results = {
        "atoms": atoms,
        "coords": ensemble["coords"],
        "energies": results,
    }
    return results


def run_gfnff_gfn2(mol, xtb_options, n_cores):
    rdDistGeom.EmbedMultipleConfs(
        mol, randomSeed=42, numConfs=100, useRandomCoords=True
    )
    if rdForceFieldHelpers.MMFFHasAllMoleculeParams(mol):
        _ = rdForceFieldHelpers.MMFFOptimizeMoleculeConfs(mol, numThreads=n_cores)
    elif rdForceFieldHelpers.UFFHasAllMoleculeParams(mol):
        _ = rdForceFieldHelpers.UFFOptimizeMoleculeConfs(mol, numThreads=n_cores)
    else:
        raise ValueError("No force field available for optimization")
    xtb_ff = ConformerCalculator(
        xtb_calculate, {"opt": None, "gfn": "ff"} | xtb_options, scr="."
    )
    xtb_gfn2_sp = ConformerCalculator(
        xtb_calculate, {"gfn": "2"} | xtb_options, scr="."
    )
    mol, _ = xtb_ff(mol, n_cores=n_cores)

    print([c.GetDoubleProp("electronic_energy") for c in mol.GetConformers()])
    # pdb.set_trace()
    mol = filter_conformers(mol, rmsdThreshold=0.5, numThreads=n_cores)
    print([c.GetDoubleProp("electronic_energy") for c in mol.GetConformers()])
    mol, _ = xtb_gfn2_sp(mol, n_cores=n_cores)
    print("here")
    print([c.GetDoubleProp("electronic_energy") for c in mol.GetConformers()])
    energies = [c.GetPropsAsDict()["electronic_energy"] for c in mol.GetConformers()]
    print(energies)
    energies = [float(c.GetProp("electronic_energy")) for c in mol.GetConformers()]
    # pdb.set_trace()
    return mol, energies


def run_gfn2_opt(atoms, coords, charge, xtb_options, n_cores):
    results = xtb_calculate(
        atoms=atoms,
        coords=coords,
        charge=charge,
        options={"gfn": "2"} | xtb_options,
        n_cores=n_cores,
    )
    return results


def generate_with_etkdg(smi, solvent, n_cores, n_confs_reactant: int = 100):
    mol3d = Chem.AddHs(Chem.MolFromSmiles(smi))

    xtb_options = {}
    if solvent:
        solvent = canonicalize_solvent(solvent, "xtb")
        xtb_options["alpb"] = solvent

    data = {"mol": mol3d, "xtb_options": xtb_options, "n_cores": n_cores}
    job = executor.submit(run_gfnff_gfn2, **data)
    gfnff_opt, energies = job.result()
    for c, e in zip(gfnff_opt.GetConformers(), energies):
        c.SetDoubleProp("electronic_energy", e)
    gfnff_opt = energy_filter_conformer(gfnff_opt, cutoff_kcalmol=15.0)
    xtb_options = {"gfn": "2", "opt": None}
    if solvent:
        xtb_options["alpb"] = canonicalize_solvent(solvent, "xtb")
    jobs = []
    atoms = [a.GetSymbol() for a in gfnff_opt.GetAtoms()]
    charge = rdmolops.GetFormalCharge(gfnff_opt)
    with executor.batch():
        for conf in gfnff_opt.GetConformers():
            coords = conf.GetPositions()
            data = {
                "atoms": atoms,
                "coords": coords,
                "charge": charge,
                "xtb_options": xtb_options,
                "n_cores": 1,
            }
            job = executor.submit(run_gfn2_opt, **data)
            jobs.append(job)
    results = [j.result() for j in jobs]
    click.echo("GFN2 Calculations completed.")
    # construct molecule
    gfn2_opt = Chem.Mol(gfnff_opt)
    gfn2_opt.RemoveAllConformers()
    for res in results:
        if res["normal_termination"]:
            conf = Chem.Conformer(len(res["coords"]))
            conf.SetPositions(np.array(res["opt_coords"]))
            conf.SetDoubleProp("electronic_energy", res["electronic_energy"])
            gfn2_opt.AddConformer(conf, assignId=True)

    gfn2_opt = energy_filter_conformer(gfn2_opt, cutoff_kcalmol=10.0)
    atoms = [a.GetSymbol() for a in gfn2_opt.GetAtoms()]
    charge = rdmolops.GetFormalCharge(gfn2_opt)
    orca_options = {"r2scan-3c": None}
    if solvent:
        orca_options["smd"] = canonicalize_solvent(solvent, "orca")
    jobs = []
    with executor.batch():
        for conf in gfn2_opt.GetConformers():
            coords = conf.GetPositions()
            data = {
                "atoms": atoms,
                "coords": coords,
                "charge": charge,
                "orca_options": orca_options,
                "n_cores": 1,
            }
            job = executor.submit(run_orca_sp, data)
            jobs.append(job)
    click.echo("Calculations submitted. Waiting for results...")
    results = [j.result() for j in jobs]
    click.echo("Calculations completed.")
    return {
        "atoms": atoms,
        "coords": [conf.GetPositions().tolist() for conf in gfn2_opt.GetConformers()],
        "energies": results,
    }


def ensemble2xyz(atoms, coords, energies):
    def ac2xyz(atoms, coords, energy):
        """Converts atom symbols and coordinates to xyz string."""
        xyz = f"{len(atoms)}\n{energy}\n"
        for atom, coord in zip(atoms, coords):
            xyz += f"{atom} {coord[0]:.8f} {coord[1]:.8f} {coord[2]:.8f}\n"
        return xyz

    xyz = ""
    for c, e in zip(coords, energies):
        xyz += ac2xyz(atoms, c, e)
        xyz += "\n"
    xyz = xyz.strip()  # Remove trailing newline
    return xyz


@click.command()
@click.option(
    "--smi",
    "-s",
    type=str,
    help="SMILES of molecule.",
)
@click.option(
    "--mode",
    type=click.Choice(["goat", "etkdg"], case_sensitive=False),
    default="goat",
    help="Calculation mode: 'goat' for GOAT calculations, 'etkdg' for ETKDG conformer generation.",
)
@click.option(
    "--n-cores",
    type=int,
    default=1,
    help="Number of CPU cores to use (default: 1).",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    default="conformers.xyz",
    help="Output file to save results (default: results.json).",
)
@click.option("--solvent", type=str, default=None, help="Solvent")
def main(smi, mode, n_cores, output, solvent):
    generate(smi, mode, n_cores, output, solvent)


def generate(smi, mode, n_cores, output, solvent):
    if mode.lower() == "goat":
        results = generate_with_goat(smi, solvent, n_cores)
    elif mode.lower() == "etkdg":
        results = generate_with_etkdg(smi, solvent, n_cores)
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'goat' or 'etkdg'.")

    # sort the coords and energies by energy
    energies = results["energies"]
    coords = results["coords"]
    atoms = results["atoms"]
    sorted_indices = np.argsort([e["electronic_energy"] for e in energies])
    results["energies"] = [energies[i] for i in sorted_indices]
    results["coords"] = [coords[i] for i in sorted_indices]
    results["atoms"] = atoms

    xyz = ensemble2xyz(**results)
    print(xyz)
    if output:
        with open(output, "w") as f:
            f.write(xyz)
    return results


if __name__ == "__main__":
    main()
    # run_gfnff_gfn2(Chem.AddHs(Chem.MolFromSmiles("c1ccccc1")), {}, 4)
