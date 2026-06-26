# AlphaFold-initialized all-atom EDG sweep code

Minimal MATLAB code for running exact and noisy all-atom EDG sweep scripts.

Before running either sweep, build the local alignment cache:

```bash
matlab -batch "run('scripts/build_mapping_cache.m')"
```

Exact sweep:

```bash
matlab -batch "run('scripts/run_exact_sweep.m')"
```

Noisy sweep:

```bash
matlab -batch "run('scripts/run_noisy_sweep.m')"
```

## Local data layout

Structural data are not committed.

Place files here:

```text
data/PDB/<PDB_ID>.pdb
data/AFDB/<AFDB_ID>.pdb
data/SIFTS/<pdb_id>.xml or <pdb_id>.xml.gz
```

The generated mapping cache is local:

```text
configs/pdb_uniprot_residue_map.csv
```

and is not committed.

## Target sheets

The exact sweep uses:

```text
configs/exact_targets.csv
```

The noisy sweep uses:

```text
configs/noisy_targets.csv
```

The state-changing examples are listed in:

```text
configs/state_targets.csv
```

To run the state-changing examples with the exact sweep, change the `targetCSV` line at the top of `scripts/run_exact_sweep.m` to:

```matlab
targetCSV = fullfile('configs', 'state_targets.csv');
```

## Outputs

Generated outputs are written under `out/` and are ignored by Git.

The exact sweep writes:

```text
out/exact_sweep/cutoff_05/K_0/<AFDB_ID>/shared/
out/exact_sweep/cutoff_05/K_0/<AFDB_ID>/runs/init_rand_seed47/
out/exact_sweep/cutoff_05/K_0/<AFDB_ID>/runs/init_randn_seed47/
out/exact_sweep/cutoff_05/K_0/<AFDB_ID>/runs/init_floyd/
out/exact_sweep/cutoff_05/K_0/<AFDB_ID>/runs/init_AF_3D/
```

The noisy sweep writes:

```text
out/noisy_sweep/cutoff_06/noise_0p003333/noise_seed_006/<AFDB_ID>/shared/
out/noisy_sweep/cutoff_06/noise_0p003333/noise_seed_006/<AFDB_ID>/runs/init_AF_rank_jitter0p001/
```

## Requirements

MATLAB with the Bioinformatics Toolbox is required because PDB parsing uses `pdbread`.
