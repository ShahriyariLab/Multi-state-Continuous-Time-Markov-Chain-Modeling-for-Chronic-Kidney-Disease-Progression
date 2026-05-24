#!/bin/bash
#SBATCH --job-name=ckd-cv-small
#SBATCH --partition=cpu
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/home/qiumengli_umass_edu/Projects/CKD/logs/ckd-cv-small-%j.out
#SBATCH --error=/home/qiumengli_umass_edu/Projects/CKD/logs/ckd-cv-small-%j.err

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
export N_REPEATS=2
export OUT_DIR=/home/qiumengli_umass_edu/Projects/CKD/results/cv_small
export MSM_FINAL_RDS=/home/qiumengli_umass_edu/Projects/CKD/results/M1_msm_hom.rds
mkdir -p "${OUT_DIR}"

Rscript scripts/ckd_msm_M1_cv.R

echo "----------------------------------------"
echo "Finished: $(date)"
