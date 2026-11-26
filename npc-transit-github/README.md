# NPC-Transit

A dependency-free Java program that estimates, from a PDB structure, the **relative passive permeability of a compact, inert protein through the nuclear
pore complex**, expressed as ln(q), the nuclear entry rate relative to a reference probe. Predictions carry an error band (leave-one-out RMSE) and a warning
for inputs outside the training range. Version 0.4.

## Status (read this first)
- **Validated, on a small dataset:** ridge regression of ln(q) on the structure-derived Stokes radius (`regress`, `predictq`). Trained on nine proteins from one study
  (Mohr et al., 2009); leave-one-out RMSE 0.71 ln-units (95% bootstrap CI 0.37–1.00) against 2.70 for a mean-only baseline. A mass-only model performs almost as well
  (0.83), so the structure-based model offers at most a small advantage. **No independent validation has been performed.**
- **Not validated:** the FG-probe scan and the coarse-grained simulation (`experimental/`). Do not use them for inference.
- **Scope:** compact, non-interacting proteins below 500 kDa; passive permeability only (no receptor-mediated transport). The tool is not designed to score viral
  capsids or other large assemblies and declines inputs above 500 kDa.

## Quick start (JDK 11+)
```
./build.sh
java -cp build npctransit.Main features my_protein.pdb
java -cp build npctransit.Main regress data/dataset_mohr2009.csv --features stokesA --out model.properties
java -cp build npctransit.Main predictq model.properties my_protein.pdb
```
Use one biological copy of the protein per file (e.g. keep chain A) and check the printed residue count and mass.

## Reproduce the manuscript results
`./reproduce.sh` builds the code, downloads the structures (`data/get_structures.sh`), runs the Stokes-radius and mass regressions, the exploratory paired comparison and
Figure 2. Outputs go to `results/`. Needs curl, awk and python3 with matplotlib.

## Layout
`src/` code · `data/` dataset and structure download script · `analysis/` post-processing scripts · `docs/MANUAL.md` manual · `experimental/` unvalidated modules ·
`test/` synthetic test data only.

## Limitations
Nine proteins from one laboratory and one assay (permeabilised HeLa cells); structure-derived radii describe deposited structures, not the labelled constructs that were assayed;
the radius estimate Rg/0.775 was below the measured Stokes radius for most proteins; wide uncertainty is expected. See the manuscript for details.

## Citation and licence
Add the paper reference and DOI here. Choose a licence before publishing (GitHub: Add file → Create new file → name it `LICENSE` → "Choose a license template").
