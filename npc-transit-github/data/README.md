# Data

`dataset_mohr2009.csv`: nine inert proteins (plus the thioredoxin reference, left unscored) from **Table 1 of Mohr et al. (2009),
EMBO J 28:2541–2553, doi:10.1038/emboj.2009.200**. Values were transcribed by hand from the published table and should be checked
against the paper. Columns: `id`, `pdb_file` (created by `get_structures.sh`), `y` = ln(q), `q` = k_i / k_thioredoxin (nuclear entry rate
constants in permeabilised HeLa cells), `q_sd` (published standard deviation of q), `stokes_nm` (measured Stokes radius), `source_doi`, `assay`, `notes`.
Lines starting with `#` are comments; a blank `y` means the entry is not scored.

PDB structures are not stored here. Run `./get_structures.sh` (needs internet). Structure-derived radii describe the deposited structures,
not the exact labelled constructs that were assayed. Cite the original paper when using these values.
