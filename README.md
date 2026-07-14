# AlphaFold-initialized all-atom EDG sweep code

This repository contains the custom MATLAB code and target configuration files used to construct matched all-atom point clouds and run the AF-EDG reconstructions reported in the manuscript.

## Repository scope

The repository contains reconstruction code, target lists, and fixed run settings. It does **not** contain third-party PDB, AlphaFold Database, or SIFTS files; generated files under `out/`; manuscript result tables; or figure source data. The numerical source values underlying the noisy plots in Fig. 2a,b are supplied with the article as **Supplementary Data 1**.

## Software environment

The manuscript analyses were run with:

- MATLAB R2025b Update 1
- Bioinformatics Toolbox 25.2
- Statistics and Machine Learning Toolbox 25.2

Bioinformatics Toolbox supplies `pdbread`, which is used for PDB parsing. Statistics and Machine Learning Toolbox supplies `procrustes`, which is used by the residue-mapping and structural-evaluation code. No other MathWorks products are required by the paper sweep scripts.

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

## Fixed settings shared by the paper sweeps

Both sweep scripts use official PDB/SIFTS-to-UniProt residue mapping, retain the full mapped residue set, construct conservatively matched all-atom point clouds with `AtomSelection="safe"`, and use `MinSeqSep=0`. PDB model 1 is used; only `ATOM` records are retained; hydrogens are removed; and blank or `A` alternate locations are accepted. They use the complete cutoff-defined graph, with no edge subsampling and no added true or AlphaFold-derived edges. Solver coordinates are saved without handedness, chirality, reflection, or other post-solver modification.

The common optimizer settings are encoded at the top of both scripts:

```matlab
opts.r           = 10000;
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

- target list: `configs/exact_targets.csv`
- all-atom cutoff: 5 Å
- unperturbed target-derived squared distances (`K=0`)
- initializations: uniform random seed 47, normal random seed 47, Floyd-Warshall, and centered AlphaFold 3D coordinates
- solver: `alternating_completion`

## State-specific reconstruction cases

The state-specific targets are listed in `configs/state_targets.csv`. Run them with the exact workflow by setting the `targetCSV` line at the top of `scripts/run_exact_sweep.m` to:

```matlab
targetCSV = fullfile('configs', 'state_targets.csv');
```

These cases use the same 5 Å, `K=0` reconstruction settings as the exact benchmark. The state-specific results reported in the manuscript are the `init_AF_3D` runs.

## Noisy reconstruction benchmark

Run:

```bash
matlab -batch "run('scripts/run_noisy_sweep.m')"
```

The noisy benchmark uses:

- target list: `configs/noisy_targets.csv`
- all-atom cutoff: 6 Å
- nominal noise levels: `eta = 0.01` and `eta = 0.05`
- implemented multiplicative scales: `K = eta/3`
- noise seeds: 3, 21, 450, 666, and 987
- initialization: AF-rank
- AF-rank embedding seed: 47
- AF-rank jitter: `1e-3`
- chemistry-aware noise: protected local chemistry edges remain unperturbed
- solver: `alternating_completion_noisy`

The script runs all 20 conditions: two proteins, two noise levels, and five seeds.

## Outputs

Both scripts write under:

```text
out/init_cutoff_sweep/
```

Exact outputs follow:

```text
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/shared/
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/runs/init_rand_seed47/
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/runs/init_randn_seed47/
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/runs/init_floyd/
out/init_cutoff_sweep/cutoff_05/<AFDB_ID>/K_0/runs/init_AF_3D/
```

Noisy outputs follow:

```text
out/init_cutoff_sweep/cutoff_06/<AFDB_ID>/K_0p00333/noise_seed_003/shared/
out/init_cutoff_sweep/cutoff_06/<AFDB_ID>/K_0p00333/noise_seed_003/runs/init_AF_rank_jitter1e-3/
```

The second noise level uses `K_0p0167`, and the remaining seed folders are `noise_seed_021`, `noise_seed_450`, `noise_seed_666`, and `noise_seed_987`.

Each script also writes:

```text
out/init_cutoff_sweep/sweep_config.mat
out/init_cutoff_sweep/manifest.csv
```

Because the exact and noisy scripts share the same output root, these two top-level summary files describe whichever sweep was run most recently. The condition-specific result folders do not overlap.

Generated outputs are ignored by Git and can be regenerated from the public structural inputs and the fixed settings above.

## Submitted code version

The manuscript should cite the full Git commit SHA corresponding to the final submitted repository state. After the last code and README commit, obtain it with:

```bash
git rev-parse HEAD
```
