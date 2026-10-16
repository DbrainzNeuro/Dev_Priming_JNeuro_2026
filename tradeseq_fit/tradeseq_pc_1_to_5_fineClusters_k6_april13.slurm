#!/bin/bash
#SBATCH --job-name=tradeseq_pcs1to5_april13_allgenes
#SBATCH --cpus-per-task=8
#SBATCH --mem=128g
#SBATCH --time=3-00:00:00
#SBATCH --gres=lscratch:150
#SBATCH --partition=norm
#SBATCH --chdir=.
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
set -x

PIPELINE_ROOT="${PIPELINE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
SCE_RDS="${SCE_RDS:-${PIPELINE_ROOT}/data/slingshot_output.rds}"

TF_GENE_LIST_RDS="${TF_GENE_LIST_RDS:-${PIPELINE_ROOT}/tradeseq_fit/gene_lists/ALL_TF.rds}"
NONTF_GENE_LIST_RDS="${NONTF_GENE_LIST_RDS:-${PIPELINE_ROOT}/tradeseq_fit/gene_lists/HVG_minus_TF.rds}"

OUTDIR="${OUTDIR:-${PIPELINE_ROOT}/results/tradeseq_fit}"
CHUNK_SIZE="${CHUNK_SIZE:-250}"
N_KNOTS="${N_KNOTS:-6}"
MIN_CELLS_EXPRESSED="${MIN_CELLS_EXPRESSED:-10}"
SEED="${SEED:-12345}"

SCRATCH_ROOT="${LSCRATCH:-/tmp}"
SCRATCH_DIR="${SCRATCH_ROOT}/tradeseq_${SLURM_JOB_ID}"
LOCAL_OUTDIR="${SCRATCH_DIR}/tradeseq_output"

module purge
module load R

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

R_SCRIPT="${PIPELINE_ROOT}/tradeseq_fit/005_run_tradeseq_biowulf.R"

mkdir -p "${SCRATCH_DIR}" "${LOCAL_OUTDIR}" "${OUTDIR}"

echo "[`date`] PIPELINE_ROOT=${PIPELINE_ROOT}"
echo "[`date`] SCE_RDS=${SCE_RDS}"
echo "[`date`] TF_GENE_LIST_RDS=${TF_GENE_LIST_RDS}"
echo "[`date`] NONTF_GENE_LIST_RDS=${NONTF_GENE_LIST_RDS}"
echo "[`date`] OUTDIR=${OUTDIR}"
echo "[`date`] SCRATCH_DIR=${SCRATCH_DIR}"
echo "[`date`] R_SCRIPT=${R_SCRIPT}"

[[ -f "${SCE_RDS}" ]] || { echo "Missing SCE_RDS: ${SCE_RDS}" >&2; exit 1; }
[[ -f "${TF_GENE_LIST_RDS}" ]] || { echo "Missing TF_GENE_LIST_RDS: ${TF_GENE_LIST_RDS}" >&2; exit 1; }
[[ -f "${NONTF_GENE_LIST_RDS}" ]] || { echo "Missing NONTF_GENE_LIST_RDS: ${NONTF_GENE_LIST_RDS}" >&2; exit 1; }
[[ -f "${R_SCRIPT}" ]] || { echo "Missing R_SCRIPT: ${R_SCRIPT}" >&2; exit 1; }

LOCAL_SCE_RDS="${SCRATCH_DIR}/$(basename "${SCE_RDS}")"
cp "${SCE_RDS}" "${LOCAL_SCE_RDS}"

LOCAL_GENE_LISTS=()
for gene_list in "${TF_GENE_LIST_RDS}" "${NONTF_GENE_LIST_RDS}"; do
  local_copy="${SCRATCH_DIR}/$(basename "${gene_list}")"
  cp "${gene_list}" "${local_copy}"
  LOCAL_GENE_LISTS+=("${local_copy}")
done

LOCAL_GENE_LIST="$(IFS=,; echo "${LOCAL_GENE_LISTS[*]}")"

Rscript "${R_SCRIPT}" \
  --sce_rds "${LOCAL_SCE_RDS}" \
  --gene_list_rds "${LOCAL_GENE_LIST}" \
  --outdir "${LOCAL_OUTDIR}" \
  --final_outdir "${OUTDIR}" \
  --chunk_size "${CHUNK_SIZE}" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --nknots "${N_KNOTS}" \
  --min_cells_expressed "${MIN_CELLS_EXPRESSED}" \
  --seed "${SEED}"

rsync -a "${LOCAL_OUTDIR}/" "${OUTDIR}/"

echo "[`date`] tradeSeq run completed successfully."
echo "[`date`] Final outputs saved in ${OUTDIR}"
