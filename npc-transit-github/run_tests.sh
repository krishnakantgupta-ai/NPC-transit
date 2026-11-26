#!/bin/sh
# Software smoke test ONLY (synthetic data; results have no scientific meaning).
java -cp build npctransit.Main train test/synthetic_dataset.csv /tmp/synth.model && \
java -cp build npctransit.Main predict /tmp/synth.model test/synth3.pdb
