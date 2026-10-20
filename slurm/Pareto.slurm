#!/usr/bin/env bash
#SBATCH --job-name=pareto_r500
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --array=0-2
#SBATCH --time=7-00:00:00
#SBATCH --output=pareto_r500_%A_%a.out
#SBATCH --error=pareto_r500_%A_%a.err

set -euo pipefail

MODEL_ROOT="/home/vinicius/ABM_HPC/mosquito"
PROJECT_ROOT="${MODEL_ROOT}/.."

cd "${MODEL_ROOT}"

strategies=(daily weekly fortnightly)
array_index="${SLURM_ARRAY_TASK_ID}"
[[ "$array_index" =~ ^[0-2]$ ]] || { echo "Invalid array index" >&2; exit 2; }
strategy="${strategies[$array_index]}"

worker_count="${PARETO_WORKERS:-31}"
if [[ ! "${worker_count}" =~ ^[0-9]+$ ]] || (( 10#${worker_count} < 1 || 10#${worker_count} > 31 )); then
    echo "PARETO_WORKERS must be an integer between 1 and 31" >&2
    exit 2
fi

export PARETO_DMW="${PARETO_DMW:-30}"
export PARETO_WASHOUT="${PARETO_WASHOUT:-30}"
export PARETO_MAX_N="${PARETO_MAX_N:-700000}"
export PARETO_RUN_SERIES="${PARETO_RUN_SERIES:-Pareto_R500}"

# One writer per strategy: protect checkpoint and CSV outputs from concurrent jobs.
lock_dir="${MODEL_ROOT}/outputs/Pareto/${strategy}"
mkdir -p "${lock_dir}"
exec 9>"${lock_dir}/.run.lock"
if ! flock -n 9; then
    echo "Another Pareto job is already writing strategy '${strategy}'." >&2
    exit 1
fi

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export GKSwstype=100

echo "Starting Pareto R=500"
echo "strategy=${strategy}"
echo "workers=${worker_count}"
echo "DMW=${PARETO_DMW} days"
echo "washout=${PARETO_WASHOUT} days"
echo "Nmax=${PARETO_MAX_N}"
echo "run series=${PARETO_RUN_SERIES}"
echo "node=${SLURMD_NODENAME:-unknown}"

srun --nodes=1 --ntasks=1 --cpus-per-task=32 --hint=nomultithread \
    /home/vinicius/software/julia-1.12.6/bin/julia \
    --project="${PROJECT_ROOT}" -p "${worker_count}" ${MODEL_ROOT}/scripts/Pareto.jl "${strategy}"
