"""Collect QM command for processing quantum mechanical calculation results."""

import json
import os

import typer
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdmolops
from tooltoad.xtb import xtb_calculate

from d2.config import NetworkConfig, IdxCounter
from d2.utils import (
    ac2xyz,
    calculate_mol_hash,
    generate_and_save,
    get_mol,
    get_random_str,
    standardize_mol,
)

app = typer.Typer()


def optimize(mol, multiplicity=1, solvent=None, n_cores=1):
    scr = os.getenv("SCRATCH", ".")
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coords = mol.GetConformer().GetPositions()
    charge = rdmolops.GetFormalCharge(mol)
    options = {"opt": None}
    if solvent is not None:
        options["alpb"] = solvent
    results = xtb_calculate(
        atoms, coords, charge, multiplicity, options=options, n_cores=n_cores, scr=scr
    )
    if not results["normal_termination"]:
        raise RuntimeError(f"XTB optimization failed:\n{results['log']}")
    opt_coords = results["opt_coords"]
    mol = Chem.MolFromXYZBlock(ac2xyz(atoms, opt_coords))
    rdDetermineBonds.DetermineBonds(mol, charge=charge)
    Chem.SanitizeMol(mol)
    for k, v in results.items():
        if "energy" in k:
            mol.SetDoubleProp(k, v)
    return [mol]

@app.command("process-md")
def process_md_command(remote: bool = False, n_cores: int = 1):
    config = NetworkConfig()
    tmp_counter = IdxCounter(str(config.network_file), "tmp_count")

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    executor = config.get_executor(slurm_job_name="md-product", cpus_per_task=n_cores)
    (config.product_data / "processed").mkdir(exist_ok=True)
    with executor.batch():
        for f in config.product_data.glob("*.json"):
            with open(f, "r") as _f:
                decoder = json.JSONDecoder()
                content = _f.read()
                data, _ = decoder.raw_decode(content)
            try:
                mol = get_mol(
                    data["atoms"],
                    data["coords"],
                    data.get("charge", 0),
                    remove_db_stereo=True,
                )
                fragments = Chem.GetMolFrags(mol, asMols=True)
                for frag in fragments:
                    Chem.SanitizeMol(frag)
                fragments = [standardize_mol(frag) for frag in fragments]
                if len(fragments) > 1:
                    sibling_hashes = [calculate_mol_hash(frag) for frag in fragments]
                    shs = "-".join([str(i) for i in sibling_hashes])
                    print(
                        f"Fragmented product {f.stem} into {len(fragments)} fragments with hashes {sibling_hashes}"
                    )
                else:
                    sibling_hashes = None
                labels = [tmp_counter.idx for _ in fragments if tmp_counter.increment()]
                for label, frag in zip(labels, fragments):

                    print(f"optimizing md product {Chem.MolToSmiles(frag)}")
                    origin_ids = [
                        int(i) for i in f.stem.split("_product")[0].split(",")
                    ]
                    properties = {"label": label}
                    if sibling_hashes:
                        properties["sibling_labels"] = [x for x in labels if x != label]
                    if remote:
                        _ = executor.submit(
                            generate_and_save,
                            optimize,
                            {
                                "mol": frag,
                                "solvent": config.solvent,
                                "n_cores": n_cores,
                            },
                            origin_ids,
                            "reaction",
                            config.new_nodes,
                            properties,
                        )

                    else:
                        generate_and_save(
                            optimize,
                            kwargs={
                                "mol": frag,
                                "solvent": config.solvent,
                                "n_cores": n_cores,
                            },
                            origin=origin_ids,
                            origin_type="reaction",
                            new_nodes_dir=config.new_nodes,
                            properties=properties,
                        )
                origin_ids.sort()
                with open(
                    config.new_reactions
                    / f"{'-'.join([str(i) for i in origin_ids])}_{'-'.join([str(label) for label in labels])}_{get_random_str()}.json",
                    "w",
                ) as new_rxn_file:
                    json.dump(
                        {
                            "reactant_ids": origin_ids,
                            "product_labels": labels,
                            "rxn_type": "md-reaction",
                        },
                        new_rxn_file,
                    )
                f.rename(config.product_data / "processed" / f.name)
                print(f"Processed {f.name} successfully.")

            except Exception as e:
                print(f"Failed to process {f}: {e}")
                continue
