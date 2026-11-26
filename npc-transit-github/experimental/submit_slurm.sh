#!/bin/bash
# Usage: slurm/submit.sh dataset.csv fg_nups.fasta   (edit partition/time for your cluster)
# One array task per (protein, epsScale) pair; results append to cache.csv. Then run `calibrate` with the same flags.
DATA=$1; FASTA=$2; CACHE=cache.csv; QUALITY=production; REPS=5; CHAINS=16
N=$(grep -v -E '^(#|id,)' "$DATA" | awk -F, '$3!=""' | wc -l)
cat > job.sbatch <<EOT
#!/bin/bash
#SBATCH --array=0-$(( N*5-1 ))
#SBATCH --cpus-per-task=$REPS
#SBATCH --time=48:00:00
#SBATCH --mem=4G
EPS=(0.5 0.75 1.0 1.5 2.0)
ROW=\$(( SLURM_ARRAY_TASK_ID / 5 )); E=\${EPS[\$(( SLURM_ARRAY_TASK_ID % 5 ))]}
LINE=\$(grep -v -E '^(#|id,)' $DATA | awk -F, '\$3!=""' | sed -n "\$(( ROW+1 ))p")
ID=\$(echo "\$LINE" | cut -d, -f1); PDB=\$(dirname $DATA)/\$(echo "\$LINE" | cut -d, -f2)
java -cp build npctransit.Main barrier "\$PDB" --id "\$ID" --eps "\$E" --fasta $FASTA --quality $QUALITY --replicas $REPS --chains $CHAINS --cache $CACHE
EOT
sbatch job.sbatch
