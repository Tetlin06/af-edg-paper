# AlphaFold-initialized all-atom EDG sweep code

This repository contains the custom MATLAB code and target configuration files used to construct matched all-atom point clouds and run the AF-EDG reconstructions reported in the manuscript.

## Repository scope

The repository contains reconstruction code, target lists, and fixed run settings. It does **not** contain third-party PDB, AlphaFold Database, or SIFTS files; generated files under `out/`; manuscript result tables; or figure source data. The numerical source values underlying the noisy plots in Fig. 2a,b are supplied with the article as **Supplementary Data 1**.

## Requirements

- MATLAB R[INSERT VERSION]
- Bioinformatics Toolbox [INSERT VERSION]
- Any additional MATLAB toolboxes actually used: [INSERT OR DELETE]

PDB parsing uses `pdbread` from the Bioinformatics Toolbox.

## Local structural-data layout

Download the public inputs from their original databases and place them as follows:

```text
data/PDB/<PDB_ID>.pdb
data/AFDB/<AFDB_ID>.pdb
data/SIFTS/<pdb_id>.xml        # or .xml.gz
```

The exact PDB IDs, AlphaFold/UniProt identifiers, and chain IDs used by each analysis are listed in:

```text
configs/exact_targets.csv
configs/noisy_targets.csv
configs/state_targets.csv
```

Build the local PDB-to-UniProt mapping cache before running a sweep:

```bash
matlab -batch "run('scripts/build_mapping_cache.m')"
```

This creates the local, untracked file:

```text
configs/pdb_uniprot_residue_map.csv
```

## Exact reconstruction sweep

```bash
matlab -batch "run('scripts/run_exact_sweep.m')"
```

## Noisy reconstruction sweep used in the manuscript

```bash
matlab -batch "run('scripts/run_noisy_sweep.m')"
```

The script runs every combination below for both rows in `configs/noisy_targets.csv`:

| Setting | Value used in the manuscript |
|---|---|
| All-atom cutoff | 6 Å |
| Nominal noise levels, eta | 0.01 and 0.05 |
| Implemented multiplicative scales, K | eta / 3 |
| Noise seeds | 3, 21, 450, 666, 987 |
| Initialization | AF-rank |
| AF-rank embedding seed | 47 |
| AF-rank jitter | 1e-3 |
| Solver rank | 10 |
| Augmented-Lagrangian parameter | 10000 |
| Maximum outer iterations | 10000 |
| Stopping tolerance | 1e-5 |
| Chemistry-aware noise | protected local chemistry edges remain unperturbed |
| BB maximum iterations | 30 |
| BB xtol | 1e-8 |
| BB gtol | 1e-8 |
| BB ftol | 1e-10 |
| BB initial alpha | 1e-3 |
| BB rho | 1e-4 |
| BB sigma | 0.1 |
| BB eta | 0.8 |

The noisy script writes a condition-specific `sweep_config.mat` and `manifest.csv` beneath each noise-level/seed folder, plus a combined file:

```text
out/noisy_sweep/manifest_all_paper_runs.csv
```

## State-specific cases

The state-specific targets are listed in `configs/state_targets.csv`. Run them with the exact reconstruction workflow using the settings reported in the manuscript.

## Outputs

Generated outputs are written under `out/` and are ignored by Git. They are reproducible from the public structural inputs and the fixed settings encoded in the scripts.

## Version used for the manuscript

At submission, record both an immutable Git tag and the full commit hash here:

```text
Commit: [INSERT FULL COMMIT HASH]
```
