# NPC-Transit 0.4: User Manual

A dependency-free Java tool that estimates, from a PDB structure, the probability that a **benign
reference protein** passes the nuclear pore complex (NPC) permeability barrier by passive diffusion.
The model is **trained by you on published experimental data**. No coefficients are shipped.

## 1. What it is and is not
- Output: a calibrated-by-validation probability P(passes), with an applicability-domain warning.
- It is a statistical descriptor model (mass, radius of gyration, hydrophobic fraction, charge density).
  It does **not** simulate FG-Nups and does **not** model receptor-mediated (importin) transport.
  An NLS-like motif is reported as an annotation only.
- A PDB structure alone cannot establish transport. All predictions are labelled "model-predicted".
- Intended scope: general/reference proteins with experimental benchmarks. Not for scoring or
  optimising viral capsids or other agents for nuclear entry.

## 2. Install / build
JDK 11+ required.
```
./build.sh                       # or: mvn package
java -cp build npctransit.Main   # prints usage
```

## 3. Commands
| Command | Purpose |
|---|---|
| `features <file.pdb>` | residues, mass, Rg, approx. Stokes radius, hydrophobic fraction, net charge, NLS-like flag |
| `train <dataset.csv> <model.properties> [lambda]` | fit logistic model; repeated stratified 5-fold CV; AUC + bootstrap CI, Brier; writes `<model>.cv_predictions.csv` |
| `predict <model.properties> <file.pdb>` | probability, domain warnings, labelled output |

Parser notes: first MODEL only, CA atoms, standard residues, altLoc A/blank. Mass = sum of residue
masses + water. Net charge = (Lys+Arg) − (Asp+Glu), pH-independent. Stokes ≈ Rg/0.775 (globular approximation).
Use the biological assembly, and remove tags/ligands consistently with how the protein was assayed.

## 4. Building the real dataset (required before any result is valid)
1. Search published NPC permeability studies (e.g. Ribbeck & Görlich 2001 EMBO J; Mohr et al. 2009 EMBO J;
   Timney et al. 2016 J Cell Biol; plus others you find). Verify each citation yourself.
2. For each protein record: PDB ID, label (1 = passes, 0 = excluded) **as measured**, DOI, assay
   (permeabilised-cell, nuclear-injection, in vitro FG-hydrogel/NPC mimic), cell system, time point, notes.
3. Use one consistent label definition (e.g. a stated fluorescence-ratio or rate threshold) and report it.
   If sources disagree, exclude the protein or document the conflict.
4. Download PDBs into `data/pdb/`, point `pdb_file` at them. Start from `data/dataset_template.csv`
   (IDs there are suggestions to verify, labels intentionally blank).
5. Expect a small dataset (tens of proteins). Keep `lambda` ≥ 0.1, report CIs, and do not claim more than the CV supports.
   The tool refuses to train with <5 examples per class.
6. Hold out an external test set (proteins from a *different* study/lab) and report performance there too.

## 5. Suggested validation for the manuscript
- Repeated stratified CV (built in) plus external hold-out; AUC with CI, Brier score, calibration plot from `cv_predictions.csv`.
- Baselines: mass-only threshold (~30–60 kDa) vs. the full model. Show whether extra features add anything.
- Ablation of each feature; sensitivity to the label threshold and lambda.
- Error analysis of misclassified proteins (e.g. highly charged, hydrophobic, or elongated/disordered).
- Limitations: small n, heterogeneous assays, no FG-Nup/receptor physics, CA-only geometry, Rg-based Stokes estimate.

## 6. Reporting text (adapt)
"NPC-Transit computes four structure-derived descriptors and fits an L2-regularised logistic regression
on N experimentally characterised proteins (sources: …). Performance was estimated by 20× repeated stratified
5-fold cross-validation and on an independent test set. Outputs are statistical estimates of passive
permeability, not mechanistic simulations or evidence of transport."

## 7. Output interpretation
`P(passes)` is only meaningful inside the training domain. If a WARNING lists features outside the training range,
do not interpret the number. Facilitated transport cannot be inferred from this tool.

## 8. Added in 0.2: docking, coarse-grained simulation, decision output
**`predict`** now prints `P(passes) = xx.x %` and a decision: YES (P >= 60%), NO (P <= 40%), INDETERMINATE (40-60%),
or NOT ASSESSED (outside training domain, or > 500 kDa; capsid-scale assemblies are out of scope by design).
The decision is only as good as the dataset you train on.

**`dock <file.pdb>`**: heuristic FG-probe docking. A Phe-ring-sized probe is scanned on a 1 A grid and scored by
hydrophobic contacts minus polar contacts, with a clash filter; top 5 clustered sites are listed with lining residues.
It is NOT a validated docking engine. Benchmark it on known FG-binding (e.g. NTF2-like, importin-beta) vs. non-binding
proteins before reporting it, and treat scores as relative only.

**`simulate <file.pdb> [epsFG] [pipeR]`**: generic coarse-grained FG-Nup brush (1 bead/residue, reduced units,
sigma = 0.6 nm, FXFG motifs, cylindrical pore, overdamped Brownian dynamics). A rigid spherical cargo (radius from the
approximate Stokes radius; F-bead attraction scaled by hydrophobic fraction) is stepped along the pore axis; mean force
is integrated to a PMF; barrier = max PMF with a 4-block SD. Prints a convergence WARNING when sampling is inadequate.
Limitations you must state: reduced pore dimensions, arbitrary epsFG and cargo mapping, short default runs, no
receptor (NTR) physics, no electrostatics. Calibrate epsFG/pipeR so that known passive cargo (e.g. GFP-sized) gives a low
barrier and large cargo a high one, run several seeds and longer trajectories, and report the calibration.
Defaults are for quick tests, not publication-grade sampling.

## 9. Not included in this package
3D visualization, UniProt/PDB web retrieval and a GUI are not part of this command-line engine, and the real
experimental reference dataset is NOT bundled: it must be curated from the primary literature (Section 4).

## 10. Version 0.3: publication workflow for the simulation
**Model.** With `--fasta fg_nups.fasta` the brush uses real FG-Nup sequences (first residue = grafted end) with the HPS
model: Ashbaugh-Hatch hydropathy interactions (epsilon0 = 0.2 kcal/mol, scaled by `epsScale`), Debye-Hueckel electrostatics
(Bjerrum 0.7 nm, Debye length 0.8 nm at ~150 mM salt), bonds 0.38 nm, sigma = 0.6 nm uniformly. Verify the lambda values in
`Hps.java` against Dignon et al. 2018 and cite the force field you use. Without `--fasta` a generic FXFG brush is used (qualitative only).
**Cargo mapping.** Radius = Rg/0.775 (Stokes approx.), lambda = mean residue lambda, charge = formal net charge. These are
approximations; state them, and ideally test alternatives (e.g. HYDROPRO Stokes radii) as a sensitivity analysis.
**Sampling.** Replicas (`--replicas`, use >= 5) are independent seeds; errors are across replicas (Student t, 95% CI).
`--quality production` = 37 windows x 1e6 steps. The tool warns when the PMF does not close.
**Calibration protocol (the part reviewers care about).**
1. Curate the dataset (Section 4) and obtain the real FG-Nup sequences (e.g. UniProt) for `fg_nups.fasta`.
2. `slurm/submit.sh dataset.csv fg_nups.fasta` runs one cluster job per (protein, epsScale); results append to `cache.csv`.
   The cache doubles as checkpointing: re-running `calibrate --cache cache.csv ...` only computes missing entries.
3. `calibrate dataset.csv --cache cache.csv --fasta ... --quality production --replicas 5 --out calibration.properties`
   reports AUC per epsScale, the selected epsScale and threshold, and **nested leave-one-out** sensitivity/specificity.
   Report the nested LOO numbers, not the in-sample ones.
4. Validate on an external set from a different laboratory/study, using the frozen `calibration.properties`
   (`simulate ... --calib calibration.properties`).
5. Sensitivity analyses: pore radius, chain number/grafting density, salt, cargo mapping, quality preset.
Keep every cache/report file and the exact command lines as supplementary material.
The test/ folder is synthetic and exists only to test the software.

## 11. Version 0.4: regression on measured permeability (recommended)
Published data (Mohr et al. 2009; Timney et al. 2016) are continuous, so version 0.4 regresses the **measured** quantity instead of
a pass/exclude label. Dataset columns: `id,pdb_file,y,q,q_sd,stokes_nm,source_doi,assay,notes`, with `y = ln(q)`,
`q = k_i / k_reference`. A blank `y` means "not scored" (use it for the reference probe, whose q = 1 by definition).
`data/dataset_mohr2009.csv` is pre-filled from Mohr Table 1 (values transcribed from the published table); you must supply/verify the
`pdb_file` entries marked in its notes. Files that do not exist are skipped with a message.

| Command | Purpose |
|---|---|
| `regress data.csv --features stokesA --lambda 0.1 --out reg.properties` | ridge regression, leave-one-out (LOO) RMSE/MAE/R2/Pearson/Spearman with bootstrap CIs, plus a mean-only baseline; also compares structure-derived vs measured Stokes radius |
| `predictq reg.properties file.pdb` | predicted ln(q) and q with +/- LOO RMSE factor and domain warning |
| `calibrate data.csv --cache c.csv --fasta fg.fasta ...` | regression mode when the header has a `y` column: scans epsScale, picks the most negative Spearman(barrier, ln q), fits ln q = a + b*barrier, reports **nested LOO** error |
| `simulate file.pdb --calib calibration.properties ...` | applies the frozen calibration (epsScale and linear map) and prints ln(q) +/- LOO RMSE |

**Rules for a defensible analysis.** Pre-specify the feature set (default: Stokes radius only; with ~9 proteins do not use more than 1-2 features).
Report LOO error relative to the mean-only baseline. Do not tune features or lambda on LOO results. Keep the thioredoxin reference out of the
scored data. Use an independent dataset (e.g. Timney et al. yeast in vivo) as the external test, converted to the same ln(q) form only if its
definition of the reference is comparable; otherwise report rank correlation (Spearman). Note the different assays (permeabilised HeLa versus
living yeast) and that the Mohr probes were chosen to be compact, inert and non-interacting, which defines the tool's domain.
