# rxnnet

A Python tool for automated exploration and visualization of chemical reaction networks using quantum mechanical calculations.

## Overview

This package enables systematic exploration of reaction networks by generating and evaluating reaction products, protomers, and tautomers through molecular dynamics simulations and DFT calculations.

## Requirements

- Python ≥3.10
- [xTB](https://github.com/grimme-lab/xtb) - for semi-empirical QM calculations
- [CREST](https://github.com/crest-lab/crest) - for tautomer/protomer generation
- [ORCA](https://orcaforum.kofo.mpg.de) - for DFT single-point energies (optional)

## Installation

```bash
pip install -e .
```

Update executable paths in your network's `config.toml` file.

## Workflow

### 1. Initialize Network

Create and configure a new reaction network:

```bash
python -m rxnnet.config <network_directory>
```

Edit the generated `config.toml` to set:
- QM parameters (temperature, pressure, solvent)
- Theory levels (xTB, DFT)
- Executable paths

### 2. Add Starting Molecule(s)

Add initial molecules from SMILES strings:

```bash
python -m rxnnet.add -s "<SMILES>"
```

### 3. Collect Molecules

Process new molecules and update the network database:

```bash
python -m rxnnet.collect
```

This assigns IDs, canonicalizes structures, and performs bookkeeping.

### 4. Generate New Structures

Expand the network using one or more generation methods:

**4a. Reaction Products (MTD simulations)**
```bash
python -m rxnnet.run_md <network_directory>/nodes/<node_id>-*.json
python -m rxnnet.process_md
```

**4b. Protomers**
```bash
python -m rxnnet.protonate <network_directory>/nodes/<node_id>-*.json -m protonate
python -m rxnnet.protonate <network_directory>/nodes/<node_id>-*.json -m deprotonate
```

**4c. Tautomers**
```bash
python -m rxnnet.tautomerize <network_directory>/nodes/<node_id>-*.json
```

### 5. Collect New Molecules

Run collect again to add newly generated structures:

```bash
python -m rxnnet.collect
```

### 6. Calculate Energies

Compute DFT single-point energies for network nodes:

```bash
python -m rxnnet.calc <network_directory>/nodes/<node_id>-*.json
```

### 7. Iterate

Repeat steps 4-6 to grow the reaction network.


### 8. Visualize

Generate an interactive HTML visualization:

```bash
python -m rxnnet.visualize
```

For follow-up analysis, see:

- [Analyzing and visualizing the reaction network](https://gist.github.com/juius/9e161982d1ae5e6cb1431e0d1bf96cf6)
- [Kinetic modelling](https://gist.github.com/juius/ac814718136c4aca7fad7f4510d40b17)

## Network Directory Structure

```
network_directory/
├── config.toml           # Configuration file
├── network-info.json     # Network metadata
├── nodes.csv             # Node database
├── nodes/                # Node data (SDF files)
├── products/             # Generated products
├── new-nodes/            # Pending nodes
├── reactions/            # Reaction data
└── qm-data/              # QM calculation results
```

## Simplified example 
```bash
# Setup
mkdir pinacol
cd pinacol
python -m rxnnet.config

# Add pinacol
python -m rxnnet.add -s "CC(C)(O)C(C)(C)O"
python -m rxnnet.collect

# Protonate
python -m rxnnet.protonate nodes/1-*.json -m protonate
python -m rxnnet.collect

# run 1 MTD simulation on protonated pinacol
python -m rxnnet.run_md nodes/2-*.json
python -m rxnnet.process_md
python -m rxnnet.collect

# run 1 MTD simulation on MTD product
python -m rxnnet.run_md nodes/3-*.json
python -m rxnnet.process_md
python -m rxnnet.collect

# Deprotonate the product of the MTD simulation
python -m rxnnet.protonate nodes/5-*.json -m deprotonate
python -m rxnnet.collect

# Calculate energies for all nodes
python -m rxnnet.calc nodes/*.json

# Filter and prune network based on variables
python -m rxnnet.status --pH 7 --temp 313.15 --max-path-energy 30 --min-count 3

# Visualize
python -m rxnnet.visualize

open network.html
```
