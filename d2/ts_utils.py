from rdkit.Chem import rdmolops
from tooltoad.chemutils import Constraint, ac2mol, get_connectivity_smiles
from tooltoad.ndscan import PotentialEnergySurface, ScanCoord
from tooltoad.orca import orca_calculate


def ts_from_scan(
    reactant,
    bond_changes: list,
    orca_options={"XTB2": None, "alpb": "water"},
    n_cores=4,
    max_cycle=10,
    verbose=True,
):
    assert len(bond_changes) == 1, "Only one bond change is allowed"
    atoms = [a.GetSymbol() for a in reactant.GetAtoms()]
    coords = reactant.GetConformer().GetPositions()
    charge = rdmolops.GetFormalCharge(reactant)
    opt_options = orca_options.copy()
    opt_options["tightopt"] = None
    preopt = orca_calculate(
        atoms,
        coords,
        charge,
        options=opt_options,
        n_cores=n_cores,
        memory=2 * n_cores,
    )
    opt_coords = preopt["opt_coords"]
    aIds = bond_changes[0][1]
    scs = [
        ScanCoord.from_current_position(
            atoms,
            opt_coords,
            aIds,
            30,
            bool(bond_changes[0][0] - 1),
        )
    ]
    scs[0].end -= 0.6
    xtb_solvent_options = {}
    if "alpb" in orca_options:
        xtb_solvent_options["alpb"] = orca_options["alpb"]
    pes = PotentialEnergySurface(atoms, opt_coords, charge, scan_coords=scs)
    pes.xtb(n_cores=n_cores, xtb_options=xtb_solvent_options, max_cycle=max_cycle)

    ts_guess = pes.traj_tensor[pes.pes_tensor.argmax()]

    if verbose:
        import shutil
        import warnings

        # Check if gnuplot is installed
        if shutil.which("gnuplot") is None:
            warnings.warn(
                "GNUplot is not installed or not in PATH. Please install it:\n"
                "  - Linux (Ubuntu/Debian): sudo apt install gnuplot\n"
                "  - macOS (Homebrew): brew install gnuplot\n"
                "  - Windows: choco install gnuplot or download from gnuplot.info",
                UserWarning,
            )
        else:
            import subprocess

            from tooltoad.chemutils import hartree2kcalmol

            # Prepare data for gnuplot
            x = pes.scan_value_tensor[:, 0]
            y = pes.pes_tensor
            y = y - y[0]
            y *= hartree2kcalmol(1)
            # Prepare the gnuplot script
            gnuplot_script = """
            set term dumb 90 30
            set xlabel "Scan Value"
            set ylabel "Electronic Energy [kcal/mol]"
            plot '-' with lines title 'PES Scan'
            """

            # Open gnuplot and send data
            process = subprocess.Popen(["gnuplot"], stdin=subprocess.PIPE, text=True)
            process.communicate(
                gnuplot_script
                + "\n".join(f"{xi} {yi}" for xi, yi in zip(x, y))
                + "\ne\n"
            )

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
            "XTB2": None,
            "alpb": "water",
            "opt": None,
        },
        xtra_inp_str=input_str,
        n_cores=n_cores,
        memory=4 * n_cores,
    )

    ts_data = ts_optimize(
        atoms,
        ts_preopt["opt_coords"],
        charge=charge,
        aIds=aIds,
        orca_options=orca_options,
        n_cores=n_cores,
    )

    ts = ac2mol(ts_data["atoms"], ts_data["opt_coords"], charge=charge)
    ts_guess = ac2mol(ts_data["atoms"], ts_data["coords"], charge=charge)
    distmat_check, max_distance = ts_distmat_check(
        ts_guess, ts, reactant, reactant, distance_threshold=0.1
    )
    freq_check = sum([d["frequency"] < 0 for d in ts_data["vibs"]]) == 1
    if verbose:
        # print number of imaginary frequencies, then the magnitude of them
        # then also if the freq check and distmat check is passed
        print(f"Imaginary frequencies: {[d['frequency'] < 0 for d in ts_data['vibs']]}")
        print(f"Max distance difference: {max_distance:.3f}")
        print(f"Distance matrix check passed: {distmat_check}")
        print(
            f"Frequency check passed: {sum([d['frequency'] < 0 for d in ts_data['vibs']]) == 1}"
        )
    return ts, ts_data, distmat_check, freq_check


def ts_optimize(
    atoms,
    coords,
    charge,
    aIds=[],
    orca_options={"XTB2": None, "alpb": "water"},
    n_cores=4,
):
    opt_options = orca_options.copy()
    opt_options["optts"] = None
    opt_options["freq"] = None
    input_str = f"""%geom
    Calc_Hess true
    TS_Mode {{ B {aIds[0]} {aIds[1]} }} end
    end"""
    ts_data = orca_calculate(
        atoms,
        coords,
        charge=charge,
        options=opt_options,
        xtra_inp_str=input_str,
        n_cores=n_cores,
        memory=4 * n_cores,
    )
    return ts_data


def ts_distmat_check(ts_guess, ts, reactant, product, distance_threshold=0.05):
    ts_dist_mat = rdmolops.Get3DDistanceMatrix(ts)
    tsguess_dist_mat = rdmolops.Get3DDistanceMatrix(ts_guess)

    adj1 = rdmolops.GetAdjacencyMatrix(reactant)
    adj2 = rdmolops.GetAdjacencyMatrix(product)
    bond_mask = adj1 | adj2
    diff_matrix = ts_dist_mat - tsguess_dist_mat
    bond_mask = bond_mask.astype(bool)
    diff_matrix[~bond_mask] = 0.0
    max_diff = abs(diff_matrix).max()
    return max_diff < distance_threshold, max_diff


def ts_irc_check(
    ts,
    reactant,
    product,
    n_cores=4,
    orca_options={"XTB2": None, "alpb": "water"},
):
    irc_options = orca_options.copy()
    irc_options["irc"] = None
    opt_options = orca_options.copy()
    opt_options["opt"] = None
    ts_atoms = [a.GetSymbol() for a in ts.GetAtoms()]
    ts_coords = ts.GetConformer().GetPositions()
    charge = rdmolops.GetFormalCharge(ts)
    irc = orca_calculate(
        ts_atoms,
        ts_coords,
        charge=charge,
        options=irc_options,
        n_cores=n_cores,
        memory=2 * n_cores,
    )
    forward = ac2mol(
        irc["irc"]["forward"]["atoms"], irc["irc"]["forward"]["opt_coords"]
    )
    backward = ac2mol(
        irc["irc"]["backward"]["atoms"], irc["irc"]["backward"]["opt_coords"]
    )

    if not get_connectivity_smiles(forward) in [
        get_connectivity_smiles(m) for m in [reactant, product]
    ]:
        print("Optimizing forward endpoint")
        forward_opt = orca_calculate(
            irc["irc"]["forward"]["atoms"],
            irc["irc"]["forward"]["opt_coords"],
            charge=charge,
            options=opt_options,
            n_cores=n_cores,
            memory=2 * n_cores,
        )
        forward = ac2mol(forward_opt["atoms"], forward_opt["opt_coords"])

    if not get_connectivity_smiles(backward) in [
        get_connectivity_smiles(m) for m in [reactant, product]
    ]:
        print("Optimizing backward endpoint")
        backward_opt = orca_calculate(
            irc["irc"]["backward"]["atoms"],
            irc["irc"]["backward"]["opt_coords"],
            charge=charge,
            options=opt_options,
            n_cores=n_cores,
            memory=2 * n_cores,
        )
        backward = ac2mol(backward_opt["atoms"], backward_opt["opt_coords"])

    irc_check = set([get_connectivity_smiles(m) for m in [forward, backward]]) == set(
        [get_connectivity_smiles(m) for m in [reactant, product]]
    )
    return irc_check, forward, backward
