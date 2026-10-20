#!/usr/bin/env bash
#SBATCH --job-name=effgain_r500
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --array=0-2
#SBATCH --time=7-00:00:00
#SBATCH --output=effgain_r500_%A_%a.out
#SBATCH --error=effgain_r500_%A_%a.err

set -euo pipefail

MODEL_ROOT="/home/vinicius/ABM_HPC/mosquito"
PROJECT_ROOT="${MODEL_ROOT}/.."

cd "${MODEL_ROOT}"

strategies=(daily weekly fortnightly)
array_index="${SLURM_ARRAY_TASK_ID}"
[[ "$array_index" =~ ^[0-2]$ ]] || {
    echo "Array index must be 0, 1, or 2" >&2
    exit 2
}
strategy="${strategies[$array_index]}"

worker_count="${EFFFREQ_WORKERS:-31}"
if [[ ! "$worker_count" =~ ^([1-9]|[12][0-9]|3[01])$ ]]; then
    echo "EFFFREQ_WORKERS must be an integer between 1 and 31" >&2
    exit 2
fi

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GKSwstype=100

output_root="${MODEL_ROOT}/outputs/EffGain"
strategy_output="${output_root}/${strategy}"
mkdir -p "${strategy_output}"
export EFFGAIN_OUTPUT_ROOT="${output_root}"

[[ -f "${MODEL_ROOT}/scripts/EffGain.jl" ]] || { echo "Missing ${MODEL_ROOT}/scripts/EffGain.jl" >&2; exit 2; }
[[ -f "${MODEL_ROOT}/scripts/main.jl" ]] || { echo "Missing ${MODEL_ROOT}/scripts/main.jl" >&2; exit 2; }

# Prevent accidental duplicate jobs for the same strategy.
exec 9>"${strategy_output}/run.lock"
if ! flock -n 9; then
    echo "Another EffGain job is already running for strategy=${strategy}" >&2
    exit 2
fi

echo "Array job: ${SLURM_ARRAY_JOB_ID:-unknown}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Strategy: ${strategy}"
echo "Node: ${SLURMD_NODENAME:-unknown}"
echo "Workers: ${worker_count}"
echo "Output: ${strategy_output}"

srun \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=32 \
    --hint=nomultithread \
    /home/vinicius/software/julia-1.12.6/bin/julia \
    --project="${PROJECT_ROOT}" \
    -p "${worker_count}" \
    "${MODEL_ROOT}/scripts/EffGain.jl" "${strategy}"
