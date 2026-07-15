# AlphaFold-initialized all-atom EDG reconstruction

This repository contains the MATLAB code and target configuration files used for the AF-EDG analyses reported in the manuscript.

## Software environment

The analyses were run with:

- MATLAB R2025b Update 1
- Bioinformatics Toolbox 25.2
- Statistics and Machine Learning Toolbox 25.2

Bioinformatics Toolbox provides `pdbread` for PDB parsing. Statistics and Machine Learning Toolbox provides `procrustes` for structural alignment and evaluation.

## Structural inputs

The scripts read experimentally determined structures from the Protein Data Bank, AlphaFold models from the AlphaFold Database, and residue mappings from SIFTS.

Arrange the files as follows:

```text
data/PDB/<PDB_ID>.pdb
data/AFDB/<AFDB_ID>.pdb
data/SIFTS/<pdb_id>.xml
```

Compressed SIFTS files with the extension `.xml.gz` are also supported.

The target identifiers and chain selections are listed in:

```text
configs/exact_targets.csv
configs/noisy_targets.csv
configs/state_targets.csv
```

Build the local PDB-to-UniProt residue-mapping cache with:

```bash
matlab -batch "run('scripts/build_mapping_cache.m')"
```

The cache is written to:

```text
configs/pdb_uniprot_residue_map.csv
```

## Shared reconstruction settings

Both reconstruction scripts use:

- SIFTS/UniProt residue correspondence
- the complete mapped residue set
- safe all-atom matching
- PDB model 1
- `ATOM` records
- blank or `A` alternate locations
- heavy atoms
- `MinSeqSep = 0`
- rank 10
- a maximum of 10,000 outer iterations
- stopping tolerance `1e-5`

The optimizer settings are:

```matlab
opts.r           = 10000;
opts.printenergy = 1;
opts.printerror  = 1;
opts.rank        = 10;
opts.maxit       = 10000;
opts.tol         = 1e-5;
opts.lamda       = opts.r;

lsopts.maxit = 30;
lsopts.xtol  = 1e-8;
lsopts.gtol  = 1e-8;
lsopts.ftol  = 1e-10;
lsopts.alpha = 1e-3;
lsopts.rho   = 1e-4;
lsopts.sigma = 0.1;
lsopts.eta   = 0.8;
```

## Exact reconstruction benchmark

Run:

```bash
matlab -batch "run('scripts/run_exact_sweep.m')"
```

The exact benchmark uses:

- `configs/exact_targets.csv`
- a 5 Å all-atom distance cutoff
- target-derived squared distances
- uniform-random initialization with seed 47
- normal-random initialization with seed 47
- Floyd-Warshall initialization
- centered AlphaFold 3D initialization
- `alternating_completion`

## State-specific reconstruction cases

The state-specific targets are listed in:

```text
configs/state_targets.csv
```

Run them with the exact reconstruction script after setting:

```matlab
targetCSV = fullfile('configs', 'state_targets.csv');
```

The manuscript reports the `init_AF_3D` reconstructions for these cases.

## Noisy reconstruction benchmark

Run:

```bash
matlab -batch "run('scripts/run_noisy_sweep.m')"
```

The noisy benchmark uses:

- `configs/noisy_targets.csv`
- a 6 Å all-atom distance cutoff
- nominal noise levels `eta = 0.01` and `eta = 0.05`
- perturbation scales `K = eta/3`
- noise seeds 3, 21, 450, 666, and 987
- chemistry-aware distance perturbations
- AF-rank initialization
- AF-rank embedding seed 47
- AF-rank jitter `1e-3`
- `alternating_completion_noisy`

The script runs 20 conditions: two proteins, two noise levels, and five seeds.

## Outputs

Both scripts write to:

```text
out/init_cutoff_sweep/
```

Exact benchmark outputs follow:

```text
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/shared/
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/runs/init_rand_seed47/
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/runs/init_randn_seed47/
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/runs/init_floyd/
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/runs/init_AF_3D/
```

Noisy benchmark outputs follow:

```text
out/init_cutoff_sweep/cutoff_06/<AFDB_ID>/K_0p00333/noise_seed_003/shared/
out/init_cutoff_sweep/cutoff_06/<AFDB_ID>/K_0p00333/noise_seed_003/runs/init_AF_rank_jitter1e-3/
```

The second noise level uses `K_0p0167`. The remaining seed directories are `noise_seed_021`, `noise_seed_450`, `noise_seed_666`, and `noise_seed_987`.

Each script also writes:

```text
out/init_cutoff_sweep/sweep_config.mat
out/init_cutoff_sweep/manifest.csv
```

The public runners save matched point clouds, retained constraints,
initial coordinates, reconstructed coordinates, and solver diagnostics.
They do not regenerate manuscript tables or figures. Reported numerical
results are provided in the manuscript tables and Supplementary Data 1.
