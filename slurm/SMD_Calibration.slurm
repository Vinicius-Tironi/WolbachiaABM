#!/usr/bin/env bash
#SBATCH --job-name=smd_calibration_r500
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
#SBATCH --output=smd_calibration_r500_%j.out
#SBATCH --error=smd_calibration_r500_%j.err

set -euo pipefail
MODEL_ROOT="/home/vinicius/ABM_HPC/mosquito"
PROJECT_ROOT="${MODEL_ROOT}/.."
SCRIPT="${MODEL_ROOT}/scripts/SMD_Calibration.jl"
[[ -f "$SCRIPT" ]] || { echo "Missing $SCRIPT" >&2; exit 2; }
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GKSwstype=100
cd "$MODEL_ROOT"
echo "Calibration=SMD_Calibration | K grid=24500:100:25500 | R=500 | 31 Julia workers | job=${SLURM_JOB_ID:-unknown}"
srun --nodes=1 --ntasks=1 --cpus-per-task=32 --hint=nomultithread \
  /home/vinicius/software/julia-1.12.6/bin/julia --project="${PROJECT_ROOT}" \
  -p 31 "$SCRIPT"
