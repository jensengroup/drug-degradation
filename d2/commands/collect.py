"""Collect command for processing new nodes and reactions."""

import shelve

import typer
from rdkit import Chem
from tooltoad.chemutils import canonicalize_resonance

from d2.config import IdxCounter, NetworkConfig
from d2.network import Node
from d2.utils import Reaction, calculate_mol_hash

app = typer.Typer()


@app.command("collect")
def collect_command():
    """Collect new nodes from new-nodes directory."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    with shelve.open(str(config.label_db), writeback=True) as db:
        if "hash2id" not in db:
            db["hash2id"] = {}
        if "label2id" not in db:
            db["label2id"] = {}
        db.sync()

        node_counter = IdxCounter(str(config.network_file), "node_count")

        # Process new SDF files
        hashes = set()
        for f in config.new_nodes.glob("*.sdf"):
            suppl = Chem.SDMolSupplier(str(f), removeHs=False)
            mol = next(suppl, None)
            mol = canonicalize_resonance(mol)
            mol_hash = calculate_mol_hash(mol)
            hashes.add(mol_hash)

            product_node_files = list(config.node_data.glob(f"*-{mol_hash}.json"))
            if len(product_node_files) == 1:
                node_idx = int(product_node_files[0].stem.split("-")[0])
                print(f"Node {mol_hash} already exists with index {node_idx}")
                node = Node.load(product_node_files[0])
                
            elif len(product_node_files) == 0:
                print(f"need to make new node with hash {mol_hash}")
                node_counter.increment()
                node = Node(
                    idx=node_counter.idx,
                    data={"gfn2-xtb": mol},
                )
                node.save(config.node_data)
                # hash2id[mol_hash] = node.idx
            else:
                print()
                print("####################################################")
                print(
                    f"Multiple nodes found for {mol_hash}, skipping: {product_node_files}"
                )
                print("####################################################")
                print()
                continue
            if mol.HasProp("label"):
                label = int(mol.GetProp("label"))
                db["label2id"][label] = node.idx
                print(f"Processing file {f} with hash {mol_hash} and label {label}")
            else:
                print(f"Processing file {f} with hash {mol_hash}")
            db["hash2id"][mol_hash] = node.idx

            db.sync()
            rxn_type = mol.GetProp("origin-type").lower()
            if rxn_type == "reaction":
                continue
            else:
                origin_idx = eval(mol.GetProp("origin-ids"))
                assert (
                    len(origin_idx) == 1
                ), "Only single origin nodes supported for now"
                origin_idx = origin_idx[0]
                # make sure reaction file doesn't exist
                rxn_file = (
                    config.reaction_data / f"{origin_idx}-{node.idx}-{rxn_type}.json"
                )
                if rxn_file.is_file():
                    print(
                        f"{rxn_type.capitalize()} reaction between {origin_idx} and {node.idx} already exists, what went wrong?"
                    )
                    continue
                else:
                    rxn = Reaction(
                        reactant_ids=[origin_idx],
                        product_ids=[node.idx],
                        rxn_type=rxn_type,
                    )
                    rxn.save(config.reaction_data)
                    print(
                        f"Created {rxn_type.capitalize()} reaction between {origin_idx} and {node.idx}"
                    )
            f.unlink()

        # Create md reactions
        for f in config.new_reactions.glob("*.json"):
            origin_ids, product_labels, _ = f.stem.split("_")
            try:
                product_ids = [
                    int(db["label2id"][int(label)])
                    for label in product_labels.split("-")
                ]
            except KeyError as e:
                print(
                    f"Error: Label {e} not found in label2id, skipping file {f}. "
                    "This might be because the product node was not created in the first pass."
                )
                continue
            product_ids.sort()
            origin_ids = [int(i) for i in origin_ids.split("-")]
            origin_ids.sort()
            if product_ids == origin_ids:
                print(f"Skipping self-reaction for {origin_ids}")
                continue
            rxn_file = (
                config.reaction_data
                / f"{'_'.join([str(i) for i in origin_ids])}-{'_'.join([str(i) for i in product_ids])}-reaction.json"
            )
            if rxn_file.is_file():
                rxn = Reaction.load(rxn_file)
                rxn.count += 1
                rxn.save(config.reaction_data)
                print(
                    f"Reaction between {origin_ids} and {product_ids} already exists, incrementing counter to {rxn.count}"
                )
            else:
                rxn = Reaction(
                    reactant_ids=origin_ids,
                    product_ids=product_ids,
                    rxn_type="reaction",
                )
                print(
                    f"Creating new reaction between {origin_ids} and {product_ids} from file {f}"
                )
                rxn.save(config.reaction_data)
            f.unlink()
