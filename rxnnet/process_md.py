"""Collect QM command for processing quantum mechanical calculation results."""

import argparse
import json
import math
import os

from rdkit import Chem
from rdkit.Chem import rdmolops
from rxnnet.xtb import xtb_calculate

from rxnnet.config import Config
from rxnnet.chemutils import (
    fragment_charge,
    get_mol,
    get_multiplicity,
    standardize_mol,
)
from rxnnet.utils import get_random_str

MAX_RECURSION_DEPTH = 5


def _preserve_atom_mapping(new_mol, old_mol):
    """Copy atom map numbers from old_mol to new_mol by atom index."""
    if new_mol is None or new_mol.GetNumAtoms() != old_mol.GetNumAtoms():
        return new_mol
    for i in range(old_mol.GetNumAtoms()):
        map_num = old_mol.GetAtomWithIdx(i).GetAtomMapNum()
        if map_num > 0:
            new_mol.GetAtomWithIdx(i).SetAtomMapNum(map_num)
    return new_mol


def _set_fragment_atom_mapping(frag_mol, parent_mol, parent_indices):
    """Set atom map numbers on fragment from parent using index mapping."""
    for frag_idx, parent_idx in enumerate(parent_indices):
        map_num = parent_mol.GetAtomWithIdx(parent_idx).GetAtomMapNum()
        if map_num > 0:
            frag_mol.GetAtomWithIdx(frag_idx).SetAtomMapNum(map_num)


def _xtb_optimize(atoms, coords, charge, multiplicity=1, solvent=None, n_cores=1):
    """Run a single XTB geometry optimization."""
    scr = os.getenv("SCRATCH", ".")
    options = {"opt": None}
    if solvent is not None:
        options["alpb"] = solvent
    results = xtb_calculate(
        atoms, coords, charge, multiplicity, options=options, n_cores=n_cores, scr=scr
    )
    if not results["normal_termination"]:
        raise RuntimeError(f"XTB optimization failed:\n{results['log']}")
    return {
        "atoms": list(atoms),
        "opt_coords": results["opt_coords"],
        **{k: v for k, v in results.items() if "energy" in k},
    }


def _determine_fragment_charges(
    fragments,
    total_charge,
    solvent=None,
    n_cores=1,
    allow_radicals=False,
):
    """Find the optimal charge (and multiplicity) assignment for a set of
    fragments."""
    scr = os.getenv("SCRATCH", ".")
    xtb_options = {}
    if solvent is not None:
        xtb_options["alpb"] = solvent

    possible_charges = fragment_charge(total_charge, len(fragments))

    if allow_radicals:
        combos_to_try = possible_charges
        print(
            f"  allow_radicals=True — evaluating all "
            f"{len(combos_to_try)} charge partition(s)."
        )
    else:
        singlet_combos = []
        radical_combos = []
        for combo in possible_charges:
            all_singlet = all(
                (sum(a.GetAtomicNum() for a in frag.GetAtoms()) - ch) % 2 == 0
                for frag, ch in zip(fragments, combo)
            )
            (singlet_combos if all_singlet else radical_combos).append(combo)

        combos_to_try = singlet_combos if singlet_combos else radical_combos
        if singlet_combos:
            print(
                f"  {len(singlet_combos)} all-singlet charge partition(s) available "
                f"(skipping {len(radical_combos)} radical partition(s))."
            )
        else:
            print(
                f"  No all-singlet partition possible; falling back to "
                f"{len(radical_combos)} radical partition(s)."
            )

    best_energy = math.inf
    best_assignment = None

    for charge_combo in combos_to_try:
        assert sum(charge_combo) == total_charge, (
            f"Charge partition {charge_combo} does not sum to {total_charge}"
        )
        system_energy = 0.0
        tmp_mults = []
        tmp_mols = []
        valid = True
        for frag, charge in zip(fragments, charge_combo):
            frag_atoms = [a.GetSymbol() for a in frag.GetAtoms()]
            frag_coords = frag.GetConformer().GetPositions()

            # Validate: can we find a valid bond ordering for this charge?
            perceived_mol = get_mol(frag_atoms, frag_coords, charge, strict=False)
            if perceived_mol is None:
                valid = False
                break
            _preserve_atom_mapping(perceived_mol, frag)
            tmp_mols.append(perceived_mol)

            multiplicity = get_multiplicity(frag, charge)
            tmp_mults.append(multiplicity)
            results = xtb_calculate(
                frag_atoms,
                frag_coords,
                charge=charge,
                multiplicity=multiplicity,
                options=xtb_options,
                n_cores=n_cores,
                scr=scr,
            )
            if results["normal_termination"]:
                system_energy += results["electronic_energy"]
            else:
                valid = False
                break
        if valid and system_energy < best_energy:
            best_energy = system_energy
            best_assignment = [
                (mol, ch, mu) for mol, ch, mu in zip(tmp_mols, charge_combo, tmp_mults)
            ]

    if best_assignment is None:
        raise RuntimeError(
            "Could not determine valid charge assignment for fragments. "
            f"Total charge {total_charge} could not be distributed over "
            f"{len(fragments)} fragment(s)."
        )

    assigned_charges = [ch for _, ch, _ in best_assignment]
    assigned_mults = [mu for _, _, mu in best_assignment]
    print(
        f"  Best partition: charges={assigned_charges}, "
        f"multiplicities={assigned_mults}, energy={best_energy:.6f} Eh"
    )
    return best_assignment


def optimize(
    mol, multiplicity=1, solvent=None, n_cores=1, allow_radicals=False, _depth=0
):
    """Optimize a molecular structure, recursively handling fragmentation."""
    if _depth > MAX_RECURSION_DEPTH:
        raise RecursionError(
            f"Recursive optimization exceeded maximum depth ({MAX_RECURSION_DEPTH}). "
            "Structure keeps fragmenting — inspect the geometry."
        )

    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coords = mol.GetConformer().GetPositions()
    charge = rdmolops.GetFormalCharge(mol)

    opt_result = _xtb_optimize(
        atoms,
        coords,
        charge,
        multiplicity=multiplicity,
        solvent=solvent,
        n_cores=n_cores,
    )

    # Perceive bonds on the optimized geometry
    opt_mol = get_mol(
        opt_result["atoms"],
        opt_result["opt_coords"],
        charge,
        strict=False,
    )
    if opt_mol is None:
        raise RuntimeError(
            f"Could not determine valid bond ordering for optimized geometry "
            f"(charge={charge}). Inspect the structure."
        )
    for k, v in opt_result.items():
        if "energy" in k:
            opt_mol.SetDoubleProp(k, v)

    _preserve_atom_mapping(opt_mol, mol)

    frag_indices = list(rdmolops.GetMolFrags(opt_mol, asMols=False))
    fragments = list(rdmolops.GetMolFrags(opt_mol, asMols=True))
    for frag_mol, atom_indices in zip(fragments, frag_indices):
        _set_fragment_atom_mapping(frag_mol, opt_mol, atom_indices)

    if len(fragments) == 1:
        return [opt_mol]

    frag_smis = [Chem.MolToSmiles(f) for f in fragments]
    print(
        f"  Optimization at depth {_depth} produced {len(fragments)} fragments "
        f"{frag_smis} (input charge={charge}, mult={multiplicity}) — "
        "determining charge partition and re-optimizing."
    )
    assignments = _determine_fragment_charges(
        fragments,
        charge,
        solvent=solvent,
        n_cores=n_cores,
        allow_radicals=allow_radicals,
    )

    # Verify total charge is conserved
    assigned_charge_sum = sum(ch for _, ch, _ in assignments)
    assert assigned_charge_sum == charge, (
        f"Charge not conserved! input={charge}, fragment sum={assigned_charge_sum}"
    )

    all_optimized = []
    for frag, frag_charge, frag_mult in assignments:
        result_frags = optimize(
            frag,
            multiplicity=frag_mult,
            solvent=solvent,
            n_cores=n_cores,
            allow_radicals=allow_radicals,
            _depth=_depth + 1,
        )
        all_optimized.extend(result_frags)

    return all_optimized


def _process_single_product(data, solvent, n_cores, allow_radicals):
    """Optimize fragments from a single MD product and return the results."""
    mol = get_mol(
        data["atoms"],
        data["opt_coords"],
        data.get("charge", 0),
        remove_db_stereo=True,
    )

    # Set initial atom map numbers (1-based for RDKit convention)
    for idx in range(mol.GetNumAtoms()):
        mol.GetAtomWithIdx(idx).SetAtomMapNum(idx + 1)

    frag_indices = list(Chem.GetMolFrags(mol, asMols=False))
    initial_fragments = list(Chem.GetMolFrags(mol, asMols=True))
    for frag_mol, atom_indices in zip(initial_fragments, frag_indices):
        _set_fragment_atom_mapping(frag_mol, mol, atom_indices)

    standardized_fragments = []
    for frag in initial_fragments:
        std_frag = standardize_mol(frag)
        _preserve_atom_mapping(std_frag, frag)
        standardized_fragments.append(std_frag)
    initial_fragments = standardized_fragments

    total_charge = data.get("charge", 0)
    total_multiplicity = data.get("multiplicity", 1)

    # Determine (fragment, multiplicity) pairs
    if len(initial_fragments) > 1:
        print(
            f"  {len(initial_fragments)} fragments "
            f"(charge={total_charge}, mult={total_multiplicity}) — "
            "determining charges before optimization."
        )
        assignments = _determine_fragment_charges(
            initial_fragments,
            total_charge,
            solvent=solvent,
            n_cores=n_cores,
            allow_radicals=allow_radicals,
        )
        frags_to_optimize = [(frag, mult) for frag, _ch, mult in assignments]
    else:
        min_mult = get_multiplicity(initial_fragments[0], total_charge)
        mult = max(total_multiplicity, min_mult)
        frags_to_optimize = [(initial_fragments[0], mult)]

    optimized = []
    for frag, mult in frags_to_optimize:
        print(f"  Optimizing {Chem.MolToSmiles(frag)} (mult={mult})")
        optimized.extend(
            optimize(
                mol=frag,
                multiplicity=mult,
                solvent=solvent,
                n_cores=n_cores,
                allow_radicals=allow_radicals,
            )
        )
    return optimized


def process_md(config: Config, n_cores: int = 1, allow_radicals: bool = False):
    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    solvent = config.settings.get("qm", {}).get("solvent")
    (config.product_data / "processed").mkdir(exist_ok=True)

    for f in config.product_data.glob("*.json"):
        with open(f, "r") as _f:
            data, _ = json.JSONDecoder().raw_decode(_f.read())
        try:
            optimized_frags = _process_single_product(
                data,
                solvent,
                n_cores,
                allow_radicals,
            )
            print(f"  {f.stem}: {len(optimized_frags)} optimized product(s)")

            products = [Chem.MolToMolBlock(frag) for frag in optimized_frags]
            mapped_product_smiles = [
                Chem.MolToSmiles(frag, canonical=False) for frag in optimized_frags
            ]
            atom_mappings = [
                [
                    atom.GetAtomMapNum() - 1
                    for atom in frag.GetAtoms()
                    if atom.GetAtomMapNum() > 0
                ]
                for frag in optimized_frags
            ]

            reaction_data = {
                "origin": data["origin"],
                "traj": data["traj"],
                "products": products,
                "mapped_products": mapped_product_smiles,
                "atom_mappings": atom_mappings,
                "reactant_atom_count": len(data["atoms"]),
            }
            with open(
                config.new_reactions
                / f"{'-'.join([str(i) for i in data['origin']])}_{get_random_str()}.json",
                "w",
            ) as new_rxn_file:
                json.dump(reaction_data, new_rxn_file)

            f.rename(config.product_data / "processed" / f.name)
            print(f"Processed {f.name} successfully.")

        except Exception as e:
            print(f"Failed to process {f}: {e}")
            continue


if __name__ == "__main__":
    """Run MD simulation on a node."""
    parser = argparse.ArgumentParser(
        description="Process reactant file into individual nodes"
    )
    parser.add_argument(
        "-d",
        "--network-dir",
        default=".",
        help="Network directory path (default: current directory)",
    )
    parser.add_argument(
        "-n", "--n-cores", type=int, default=1, help="Number of cores (default: 1)"
    )
    parser.add_argument(
        "--allow-radicals",
        action="store_true",
        default=False,
        help="Consider radical (open-shell) charge partitions when fragments form",
    )

    args = parser.parse_args()
    config = Config(args.network_dir)
    process_md(config, n_cores=args.n_cores, allow_radicals=args.allow_radicals)
