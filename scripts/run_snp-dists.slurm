#!/bin/bash
#SBATCH --job-name=snpdists
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# Usage: sbatch submit.slurm PREFIX [ROARY_DIR] [OUTDIR]
if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch submit.slurm PREFIX [ROARY_DIR] [OUTDIR]"
  exit 1
fi

PREFIX="$1"
ROARY_DIR="${2:-$PWD}"
OUTDIR="${3:-$ROARY_DIR}"

mkdir -p "$OUTDIR"

# Canonicalize paths
ALN="$ROARY_DIR/core_gene_alignment.aln"
if [[ ! -s "$ALN" ]]; then
  echo "ERROR: Alignment not found or empty: $ALN"
  exit 2
fi

echo "Job started: $(date)"
echo "Prefix:      $PREFIX"
echo "ROARY_DIR:   $(readlink -f "$ROARY_DIR" 2>/dev/null || echo "$ROARY_DIR")"
echo "OUTDIR:      $(readlink -f "$OUTDIR" 2>/dev/null || echo "$OUTDIR")"
echo "Alignment:   $ALN"

# --- Env activation (tries mamba, then conda, else modules) ---
if command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook --shell=bash)"
  mamba activate snp-dists || true
elif command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate snp-dists || true
fi
command -v snp-dists >/dev/null || { module load snp-dists 2>/dev/null || true; }

echo "snp-dists path: $(command -v snp-dists || echo 'NOT FOUND')"
snp-dists -v

set -x  # log commands
# Wide matrix (TSV) with header
snp-dists -q -b "$ALN" > "${OUTDIR}/snpdist_${PREFIX}.tsv"
# Long (molten) 3-column format
snp-dists -q -m "$ALN" > "${OUTDIR}/snpdist_${PREFIX}_long.tsv"
set +x

# Report sizes so the log proves files exist
echo "Outputs:"
ls -lh "${OUTDIR}/snpdist_${PREFIX}.tsv" "${OUTDIR}/snpdist_${PREFIX}_long.tsv"
echo "First lines:"
head -n 3 "${OUTDIR}/snpdist_${PREFIX}.tsv"

echo "Job ended: $(date)"

