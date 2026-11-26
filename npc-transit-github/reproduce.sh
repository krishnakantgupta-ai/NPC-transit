#!/bin/sh
# Reproduces the regression results of the manuscript. Needs: JDK 11+, curl, awk, python3 with matplotlib.
set -e
cd "$(dirname "$0")"
./build.sh
./data/get_structures.sh
mkdir -p results
J="java -cp build npctransit.Main"
$J regress data/dataset_mohr2009.csv --features stokesA --out results/reg_stokes.properties | tee results/run_stokes.txt
$J regress data/dataset_mohr2009.csv --features massKDa --out results/reg_mass.properties | tee results/run_mass.txt
python3 analysis/paired_comparison.py results/reg_stokes.properties.loo_predictions.csv results/reg_mass.properties.loo_predictions.csv | tee results/paired_comparison.txt
python3 analysis/make_figure2.py results/reg_stokes.properties.loo_predictions.csv results/figure2.png
