#!/bin/bash
#SBATCH --job-name=gwas_glm
#SBATCH --output=logs/gwas_glm.%A.%a.out
#SBATCH --error=logs/gwas_glm.%A.%a.err
#SBATCH -p kellis
#SBATCH --array=1-50
#SBATCH --time=2-00:00:00
#SBATCH --mem=12000
#SBATCH --cpus-per-task=2

set -euo pipefail

job_list="glm_job_lists/all_jobs.sh"
total_jobs=$(wc -l < "$job_list")
total_arrays=50

base_batch=$(( total_jobs / total_arrays ))       # 11
remainder=$(( total_jobs % total_arrays ))         # 2

if [ "$SLURM_ARRAY_TASK_ID" -le "$remainder" ]; then
    batch_size=$(( base_batch + 1 ))
    start_line=$(( (SLURM_ARRAY_TASK_ID - 1) * batch_size + 1 ))
else
    batch_size=$base_batch
    start_line=$(( remainder * (base_batch + 1) + (SLURM_ARRAY_TASK_ID - remainder - 1) * base_batch + 1 ))
fi

end_line=$(( start_line + batch_size - 1 ))

# Guard against overrun
if [ "$start_line" -gt "$total_jobs" ]; then
    echo "Task $SLURM_ARRAY_TASK_ID has no lines to run."
    exit 0
fi
if [ "$end_line" -gt "$total_jobs" ]; then
    end_line=$total_jobs
fi

echo "Array task $SLURM_ARRAY_TASK_ID running lines $start_line to $end_line" >> logs/array_batch_ranges.log

for i in $(seq $start_line $end_line); do
  CMD=$(sed -n "${i}p" "$job_list")
  echo "[$SLURM_ARRAY_TASK_ID:$i] Running: $CMD"
  eval "$CMD"
done

