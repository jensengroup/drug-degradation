from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdmolops
from tooltoad.orca import orca_calculate
from d2.network import *
import logging
import os
from tooltoad.md import md_step
from tooltoad.xtb import xtb_calculate


def generate_initial_conformer(smi: str):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"Could not interpret input data as SMILES: {smi}")
    mol3d = Chem.AddHs(mol)
    _ = rdDistGeom.EmbedMolecule(mol3d)
    atoms = [a.GetSymbol() for a in mol3d.GetAtoms()]
    coords = mol3d.GetConformer().GetPositions()
    charge = rdmolops.GetFormalCharge(mol3d)
    return atoms, coords, charge


def generate_best_conformer(atoms, coords, charge, mutliplicity, options, scr):
    pass


def calculate_node(
    atoms, coords, charge, multiplicity, options, lot, scr, n_cores, idx, node_data_dir
):
    results = orca_calculate(
        atoms=atoms,
        coords=coords,
        charge=charge,
        multiplicity=multiplicity,
        options=options,
        scr=scr,
        n_cores=n_cores,
    )
    node = Node(type="species", data=[QMResult({lot: results})], idx=idx)
    node.save(node_data_dir / f"{node.idx}-{node.id}.json")


def run_md(
    atoms,
    coords,
    charge,
    multiplicity,
    options,
    detailed_str,
    n_cores,
    product_data_dir,
    species_ids,
    unique_id,
):
    _logger = logging.getLogger("tooltoad.md")
    _logger.setLevel(logging.INFO)
    _logger.addHandler(logging.StreamHandler())
    SCR = os.getenv("SCRATCH", ".")
    products = md_step(
        atoms,
        coords,
        charge,
        multiplicity,
        options=options,
        detailed_input_str=detailed_str,
        n_md_cores=n_cores,
        max_products=1,
        scr=SCR,
    )
    # only keep first product, TODO: potentially keep first N products later on
    if len(products) > 1:
        products = [products[0]]
    print("MD finished")
    for frame, product in products:
        product = orca_calculate(
            product["atoms"],
            product["coords"],
            charge=charge,
            multiplicity=multiplicity,
            scr=SCR,
            options=options,
            n_cores=n_cores,
        )
        product["frame"] = frame
        with open(
            product_data_dir / f"{species_ids}_{unique_id}_product.json", "w"
        ) as f:
            json.dump(product, f)
