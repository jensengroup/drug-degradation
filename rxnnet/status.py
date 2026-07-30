"""Filtered network tables for the CLI.

The node/edge filtering here mirrors Data.parseData()/Data.filter() in
templates/js/data.js exactly (barrier propagation from a hyperedge's primary
product to its eliminated byproducts, junction/hyperedge pruning, edge count
threshold) so that `--max-path-energy`/`--min-count` produce the same
"Nodes: N" / "Edges: N" counts as the matching sliders in the interactive
visualization.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from rdkit import Chem

from rxnnet.network import ReactionNetwork


def _clean_smiles(mol: Optional[Chem.Mol]) -> str:
    """Canonical SMILES with explicit Hs and atom maps stripped, for display."""
    if mol is None:
        return ""
    m = Chem.RemoveHs(mol)
    for atom in m.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(m)


def _propagate_eliminated_barriers(
    nodes: Dict[int, Dict[str, Any]], original_edges: List[Dict[str, Any]]
) -> None:
    """Give eliminated byproducts the barrier of their hyperedge's primary product.

    Eliminated products never get their own pathway traced by
    compute_pathways() (it only follows primary products), so they'd
    otherwise always pass the filter regardless of the cutoff.
    """
    for e in original_edges:
        eliminated = e.get("smaller_products") or []
        if not eliminated:
            continue
        primary = nodes.get(e["end"])
        if primary is None:
            continue
        for elim_id in eliminated:
            elim_node = nodes.get(elim_id)
            if elim_node is None:
                continue
            if elim_node["path_energy"] is None or (
                primary["path_energy"] is not None
                and primary["path_energy"] < elim_node["path_energy"]
            ):
                elim_node["path_energy"] = primary["path_energy"]


def _expand_edges(original_edges: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Split hyperedges into reactant->junction->product legs, same as the JS Data.parseData()."""
    edges = []
    junction_id = -1
    for e in original_edges:
        eliminated = e.get("smaller_products") or []
        count = e.get("count", 1)
        rxn_type = e.get("type", "")

        if eliminated:
            jid = junction_id
            junction_id -= 1
            edges.append(
                {"from": e["begin"], "to": jid, "type": rxn_type, "count": count, "kind": "to_junction"}
            )
            edges.append(
                {"from": jid, "to": e["end"], "type": rxn_type, "count": count, "kind": "primary"}
            )
            for elim_id in eliminated:
                edges.append(
                    {"from": jid, "to": elim_id, "type": rxn_type, "count": count, "kind": "eliminated"}
                )
        else:
            edges.append(
                {"from": e["begin"], "to": e["end"], "type": rxn_type, "count": count, "kind": "simple"}
            )
    return edges


def filter_network_data(
    data: Dict[str, Any],
    max_path_energy: float = 50.0,
    min_count: int = 1,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Apply the same filtering as the visualization's sliders to raw visualization data.

    `data` is the dict returned by ReactionNetwork.to_visualization_data().
    Returns (filtered_nodes, filtered_edges) as lists of plain dicts.
    """
    substrate_id = data["substrateId"]
    rel_energy_map = data["relEnergyMap"]
    charge_map = data["chargeMap"]
    weight_map = data["weightMap"]
    best_barrier = data["bestPathBarrierMap"]
    original_edges = data["originalEdgeData"]

    nodes: Dict[int, Dict[str, Any]] = {}
    for nid in data["svgMap"].keys():
        nodes[nid] = {
            "id": nid,
            "is_substrate": nid == substrate_id,
            "rel_energy": rel_energy_map.get(nid, 0.0),
            "path_energy": best_barrier.get(nid),
            "charge": charge_map.get(nid, 0),
            "weight": weight_map.get(nid, 0),
        }

    _propagate_eliminated_barriers(nodes, original_edges)

    visible_node_ids = {
        nid
        for nid, n in nodes.items()
        if n["is_substrate"]
        or n["path_energy"] is None
        or n["path_energy"] <= max_path_energy
    }

    edges = _expand_edges(original_edges)

    def edge_passes(e: Dict[str, Any]) -> bool:
        count_pass = (e["count"] or 1) >= min_count
        if e["kind"] == "to_junction":
            return count_pass and e["from"] in visible_node_ids
        if e["kind"] in ("primary", "eliminated"):
            return count_pass and e["to"] in visible_node_ids
        return count_pass and e["from"] in visible_node_ids and e["to"] in visible_node_ids

    filtered_edges = [e for e in edges if edge_passes(e)]

    # A hyperedge junction only survives if it still has a reactant coming in
    # and at least one product going out; drop dangling stub legs otherwise.
    junction_incoming = {e["to"] for e in filtered_edges if e["kind"] == "to_junction"}
    junction_outgoing = {
        e["from"] for e in filtered_edges if e["kind"] in ("primary", "eliminated")
    }
    valid_junctions = junction_incoming & junction_outgoing

    filtered_edges = [
        e
        for e in filtered_edges
        if not (e["kind"] == "to_junction" and e["to"] not in valid_junctions)
        and not (e["kind"] in ("primary", "eliminated") and e["from"] not in valid_junctions)
    ]

    connected_ids = set()
    for e in filtered_edges:
        if e["from"] > 0:
            connected_ids.add(e["from"])
        if e["to"] > 0:
            connected_ids.add(e["to"])
    if substrate_id:
        connected_ids.add(substrate_id)

    filtered_nodes = [nodes[nid] for nid in nodes if nid in connected_ids]
    filtered_nodes.sort(key=lambda n: n["id"])

    # Junctions are just a rendering artifact of the graph (the split point
    # for a hyperedge's multiple products). For the table, collapse each leg
    # back into a direct reactant -> product edge instead of exposing the
    # synthetic junction id.
    junction_reactant = {
        e["to"]: e["from"] for e in filtered_edges if e["kind"] == "to_junction"
    }
    display_edges = []
    for e in filtered_edges:
        if e["kind"] == "to_junction":
            continue
        if e["kind"] in ("primary", "eliminated"):
            display_edges.append({**e, "from": junction_reactant[e["from"]]})
        else:
            display_edges.append(e)
    display_edges.sort(key=lambda e: (e["from"], e["to"]))

    return filtered_nodes, display_edges


def build_tables(
    network_dir: str | Path,
    substrate_id: int = 1,
    pH: float = 7.0,
    temperature: float = 313.15,
    max_path_energy: float = 50.0,
    min_count: int = 1,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load a network, filter it, and return (nodes_df, edges_df)."""
    network = ReactionNetwork(
        network_dir=network_dir,
        substrate_id=substrate_id,
        pH=pH,
        temperature=temperature,
    )
    network.compute_pathways()

    data = network.to_visualization_data()
    nodes, edges = filter_network_data(
        data, max_path_energy=max_path_energy, min_count=min_count
    )
    for n in nodes:
        mol_node = network.nodes.get(n["id"])
        n["smiles"] = _clean_smiles(mol_node.mol) if mol_node else ""

    nodes_df = pd.DataFrame(
        nodes,
        columns=["id", "rel_energy", "path_energy", "charge", "weight", "smiles"],
    )
    edges_df = pd.DataFrame(edges, columns=["from", "to", "type", "count"])

    return nodes_df, edges_df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Print a filtered view of the reaction network as a table, "
        "using the same filtering as the interactive visualization."
    )
    parser.add_argument(
        "-d",
        "--network-dir",
        default=".",
        help="Network directory path (default: current directory)",
    )
    parser.add_argument("-s", "--substrate", default=1, type=int, help="Substrate ID")
    parser.add_argument("--pH", type=float, default=7.0, help="pH")
    parser.add_argument("--temp", type=float, default=313.15, help="Temperature (K)")
    parser.add_argument(
        "--max-path-energy",
        type=float,
        default=50.0,
        help="Max path energy (barrier) cutoff in kcal/mol, same as the "
        "visualization's 'Max Path Energy' slider (default: 50)",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Minimum reaction count threshold, same as the visualization's "
        "'Edge Count Threshold' slider (default: 1)",
    )
    parser.add_argument(
        "--mode",
        choices=["nodes", "edges", "both"],
        default="both",
        help="Which table(s) to print (default: both)",
    )
    parser.add_argument(
        "--csv",
        help="Write CSV instead of printing. With --mode both, used as a "
        "prefix: <csv>-nodes.csv and <csv>-edges.csv",
    )

    args = parser.parse_args()

    nodes_df, edges_df = build_tables(
        network_dir=args.network_dir,
        substrate_id=args.substrate,
        pH=args.pH,
        temperature=args.temp,
        max_path_energy=args.max_path_energy,
        min_count=args.min_count,
    )

    def emit(name: str, df: pd.DataFrame) -> None:
        if args.csv:
            path = (
                Path(args.csv)
                if args.mode != "both"
                else Path(f"{args.csv}-{name}.csv")
            )
            df.to_csv(path, index=False)
            print(f"Wrote {len(df)} {name} to {path}")
        else:
            print(f"\n{name.upper()} ({len(df)})")
            print(df.to_string(index=False))

    if args.mode in ("nodes", "both"):
        emit("nodes", nodes_df)
    if args.mode in ("edges", "both"):
        emit("edges", edges_df)
