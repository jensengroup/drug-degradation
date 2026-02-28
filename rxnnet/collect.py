"""Collect command for processing new nodes and reactions."""

import argparse
import json

import pandas as pd
from rdkit import Chem
from rxnnet.chemutils import canonicalize_resonance, get_smiles

from rxnnet.config import Config, IdxCounter
from rxnnet.node import Node
import logging

_logger = logging.getLogger(__name__)


def _load_node_mol(node_data_dir, node_id):
    """Load RDKit mol from node JSON file."""
    matches = sorted(node_data_dir.glob(f"{node_id}-*.json"))
    if not matches:
        raise FileNotFoundError(f"No node file found for id {node_id}")
    node = Node.load(matches[0])
    return next(iter(node._data.values()))


def _assign_sequential_atom_maps(mols):
    """Assign sequential atom map numbers across multiple molecules.

    Returns list of SMILES with atom mapping.
    """
    global_map = 1
    mapped_smiles = []
    for mol in mols:
        mol_copy = Chem.Mol(mol)
        for atom in mol_copy.GetAtoms():
            atom.SetAtomMapNum(global_map)
            global_map += 1
        mapped_smiles.append(Chem.MolToSmiles(mol_copy, canonical=False))
    return mapped_smiles


def collect(config):
    """Collect new nodes from new-nodes directory."""

    # load nodes.csv
    nodes_file = config.network_dir / "nodes.csv"
    if not nodes_file.exists():
        nodes_df = pd.DataFrame(columns=["idx", "canonical_smiles"])
        nodes_df.to_csv(nodes_file, index=False)

    nodes_df = pd.read_csv(nodes_file)
    existing_nodes = {
        smi: int(idx)
        for _, (idx, smi) in nodes_df[["idx", "canonical_smiles"]].iterrows()
    }
    node_counter = IdxCounter(str(config.network_file), "node_count")

    # load reactions.csv
    reactions_file = config.network_dir / "reactions.csv"
    if not reactions_file.exists():
        reactions_df = pd.DataFrame(
            columns=[
                "reactant_ids",
                "product_ids",
                "reaction_smiles",
                "mapped_reaction_smiles",
                "count",
                "rxn_type",
            ]
        )
        reactions_df.to_csv(reactions_file, index=False)

    reactions_df = pd.read_csv(reactions_file)
    if "mapped_reaction_smiles" not in reactions_df.columns:
        reactions_df["mapped_reaction_smiles"] = reactions_df.get("reaction_smiles", "")
    if "reaction_smiles" not in reactions_df.columns:
        reactions_df["reaction_smiles"] = reactions_df.get("mapped_reaction_smiles", "")

    existing_reactions = {}
    for _, row in reactions_df.iterrows():
        mapped_reaction_smiles = row.get("mapped_reaction_smiles", "")
        existing_reactions[mapped_reaction_smiles] = {
            "reactant_ids": row.get("reactant_ids"),
            "product_ids": row.get("product_ids"),
            "reaction_smiles": row.get("reaction_smiles", mapped_reaction_smiles),
            "count": row.get("count", 1),
            "rxn_type": row.get("rxn_type", ""),
        }

    processed_nodes = []
    processed_reactions = []
    processed_nodes_dir = config.new_nodes / "processed"
    processed_reactions_dir = config.new_reactions / "processed"
    for d in [processed_nodes_dir, processed_reactions_dir]:
        d.mkdir(exist_ok=True)

    # for protonation/tautomer noders
    for f in config.new_nodes.glob("*.sdf"):
        suppl = Chem.SDMolSupplier(str(f), removeHs=False)
        mol = next(suppl, None)
        mol = canonicalize_resonance(mol)

        smiles = get_smiles(mol)
        if smiles not in existing_nodes:
            print("new smiles")
            # write a node file, this can then store 3d information
            node = Node(
                idx=node_counter.idx,
                data={"gfn2-xtb": mol},
            )
            node.save(config.node_data)
            node_idx = node_counter.idx
            existing_nodes[smiles] = int(node_idx)
            _logger.info(f"Added new node {node_idx} for SMILES {smiles}")
            node_counter.increment()
        # update reactions.csv
        idx2smiles = {v: k for k, v in existing_nodes.items()}
        reactant_idx = eval(mol.GetProp("origin-ids"))[0]
        rxn_type = mol.GetProp("origin-type").lower()
        # get smiles for each product
        product_smiles = get_smiles(mol)
        # get smiles for reactants
        reactant_smiles = idx2smiles[reactant_idx] if rxn_type != "add" else ""
        # make reaction smiles
        rxn_smiles = f"{reactant_smiles}>>{product_smiles}"
        mapped_rxn_smiles = rxn_smiles
        if mapped_rxn_smiles in existing_reactions:
            print(rxn_smiles)
            print()
            print(existing_reactions)
            # import pdb
            # pdb.set_trace()
        else:
            # make new reaction entry
            existing_reactions[mapped_rxn_smiles] = {
                "reactant_ids": [reactant_idx],
                "product_ids": [existing_nodes[smiles]],
                "reaction_smiles": rxn_smiles,
                "count": 1,
                "rxn_type": rxn_type,
            }
            _logger.info(f"Added new reaction {rxn_smiles}")
        processed_nodes.append(f)

    # for mtd products
    idx2smiles = {v: k for k, v in existing_nodes.items()}
    for f in config.new_reactions.glob("*.json"):
        print("processing", f)
        with open(f, "r") as _f:
            data = json.load(_f)
        reactant_ids = data["origin"]
        products = [
            Chem.MolFromMolBlock(block, removeHs=False) for block in data["products"]
        ]
        mapped_products = data.get("mapped_products", [])

        # add new products
        product_ids = []
        for product in products:
            smiles = get_smiles(product)
            if smiles not in existing_nodes:
                node = Node(
                    idx=node_counter.idx,
                    data={"gfn2-xtb": product},
                )
                node.save(config.node_data)
                node_idx = node_counter.idx
                existing_nodes[smiles] = int(node_idx)
                product_ids.append(int(node_idx))
                node_counter.increment()
                _logger.info(f"Added new node {node_idx} for SMILES {smiles}")
            else:
                product_ids.append(existing_nodes[smiles])
        # get smiles for each product
        product_smiles = [get_smiles(mol) for mol in products]
        product_smiles.sort(key=lambda x: len(x), reverse=True)
        # get smiles for reactants
        reactant_smiles = [idx2smiles[idx] for idx in reactant_ids]
        reactant_smiles.sort(key=lambda x: len(x), reverse=True)
        # make reaction smiles
        rxn_smiles = f"{'.'.join(reactant_smiles)}>>{'.'.join(product_smiles)}"

        # Generate mapped reaction SMILES
        mapped_rxn_smiles = rxn_smiles  # Default to unmapped
        if mapped_products and len(mapped_products) == len(products):
            # Use pre-generated mapped product SMILES from process_md
            reactant_mols = [
                _load_node_mol(config.node_data, idx) for idx in reactant_ids
            ]
            mapped_reactant_smiles = _assign_sequential_atom_maps(reactant_mols)
            # Sort to match canonical ordering
            mapped_products_sorted = sorted(
                zip(mapped_products, product_smiles),
                key=lambda x: len(x[1]),
                reverse=True,
            )
            mapped_product_smiles = [m for m, _ in mapped_products_sorted]
            mapped_rxn_smiles = (
                f"{'.'.join(mapped_reactant_smiles)}>>{'.'.join(mapped_product_smiles)}"
            )

        if mapped_rxn_smiles in existing_reactions:
            # update count
            existing_reactions[mapped_rxn_smiles]["count"] += 1
            _logger.info(
                "Updated reaction count to "
                f"{existing_reactions[mapped_rxn_smiles]['count']} for {rxn_smiles}"
            )
        else:
            # make new reaction entry
            existing_reactions[mapped_rxn_smiles] = {
                "reactant_ids": reactant_ids,
                "product_ids": product_ids,
                "reaction_smiles": rxn_smiles,
                "count": 1,
                "rxn_type": "mtd-reaction",
            }
            _logger.info(f"Added new reaction {rxn_smiles}")
        processed_reactions.append(f)
    # update nodes.csv and reactions.csv
    nodes_df = pd.DataFrame(
        {
            "idx": list(existing_nodes.values()),
            "canonical_smiles": list(existing_nodes.keys()),
        }
    )
    nodes_df.to_csv(nodes_file, index=False)

    reactions_df = pd.DataFrame(
        {
            "reactant_ids": [v["reactant_ids"] for v in existing_reactions.values()],
            "product_ids": [v["product_ids"] for v in existing_reactions.values()],
            "reaction_smiles": [
                v["reaction_smiles"] for v in existing_reactions.values()
            ],
            "mapped_reaction_smiles": list(existing_reactions.keys()),
            "count": [v["count"] for v in existing_reactions.values()],
            "rxn_type": [v["rxn_type"] for v in existing_reactions.values()],
        }
    )
    reactions_df.to_csv(reactions_file, index=False)

    # move processed files
    for f in processed_nodes:
        f.rename(processed_nodes_dir / f.name)
    for f in processed_reactions:
        f.rename(processed_reactions_dir / f.name)


def main():
    """Collect."""
    parser = argparse.ArgumentParser(description="Collect")
    parser.add_argument(
        "-d",
        "--network-dir",
        default=".",
        help="Network directory path (default: current directory)",
    )

    args = parser.parse_args()

    config = Config(args.network_dir)

    if not config.is_initialized():
        print("Network not initialized. Run setup first.")
        return

    collect(config)


if __name__ == "__main__":
    logger = logging.getLogger("collect")
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler())
    main()
