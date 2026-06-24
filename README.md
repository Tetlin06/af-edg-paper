# AlphaFold-initialized all-atom EDG paper code

Minimal MATLAB code for reproducing the paper experiments.

This repo intentionally excludes exploratory code: no subsampling, no clustering, no added edges, no PAE/pLDDT edge selection, no PyMOL export, and no exploratory main script.

Structural data are not committed.

Expected local data layout:

data/PDB/<PDB_ID>.pdb
data/AFDB/<AFDB_ID>.pdb
data/SIFTS/<pdb_id>.xml or <pdb_id>.xml.gz

Build alignment cache:

matlab -batch "run('scripts/build_mapping_cache.m')"

Paper entry points:

matlab -batch "run('scripts/reproduce_exact_constraints.m')"
matlab -batch "run('scripts/reproduce_noisy_constraints.m')"
matlab -batch "run('scripts/reproduce_state_cases.m')"
