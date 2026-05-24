#!/bin/bash
#SBATCH --job-name=ckd-cv
#SBATCH --partition=cpu
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=192G
#SBATCH --time=02:00:00
#SBATCH --constraint=amd7763
#SBATCH --output=/home/qiumengli_umass_edu/Projects/CKD/logs/ckd-cv-%j.out
#SBATCH --error=/home/qiumengli_umass_edu/Projects/CKD/logs/ckd-cv-%j.err

set -euo pipefail

source /usr/share/lmod/lmod/init/bash
module load conda/latest
conda activate /work/pi_lshahriyari_umass_edu/qiumengli_umass_edu/ckd_env

cd /home/qiumengli_umass_edu/Projects/CKD

echo "Host:    $(hostname)"
echo "CPUs:    ${SLURM_CPUS_PER_TASK}"
echo "Started: $(date)"
echo "----------------------------------------"

export MC_CORES="${SLURM_CPUS_PER_TASK}"
export N_REPEATS=10

# Prevent BLAS thread oversubscription under mclapply fork parallelism.
# Without this, each of MC_CORES workers spawns ~128 BLAS threads,
# crushing the node and stalling every msm fit.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

Rscript scripts/ckd_msm_M1_cv.R

echo "----------------------------------------"
echo "Finished: $(date)"
