#!/bin/bash
#SBATCH --job-name=roary_fasttree
#SBATCH --partition=cpu2023
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/ST131_E.coli/fasta_links
#SBATCH --output=/work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/ST131_E.coli/fasta_links/gff_files/logs/%x_%j.out
#SBATCH --error=/work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/ST131_E.coli/fasta_links/gff_files/logs/%x_%j.err
#SBATCH --open-mode=append

set -Eeuo pipefail
set -x

# ---- paths ----
PROJECT_DIR="/work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/ST131_E.coli/fasta_links"
GFF_DIR="${PROJECT_DIR}/gff_files"
OUTDIR_BASE="${PROJECT_DIR}/gff_files/roary_out"    # <-- do NOT mkdir this; Roary will create
THREADS="${SLURM_CPUS_PER_TASK:-32}"
CORE_PCT=95
ID_PCT=95

mkdir -p "${GFF_DIR}/logs"

# ---- environment bootstrap ----
module purge || true
export MAMBA_ROOT_PREFIX=/work/vetmed_data/mamba
eval "$(mamba shell hook --shell bash)"
mamba activate roary_fasttree

# Clean Perl + locale
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
export PERL5LIB=""
export LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8

# Sanity checks (fixed: -w shows version; -v = verbose)
perl -MFile::Find::Rule -e1 && echo "File::Find::Rule OK"
roary -w

# --- gather inputs robustly ---
shopt -s nullglob
GFFS=( "${GFF_DIR}"/*.gff )
echo "Found ${#GFFS[@]} GFF files under ${GFF_DIR}"
if (( ${#GFFS[@]} < 2 )); then
  echo "ERROR: need ≥2 .gff files; got ${#GFFS[@]}"; exit 2
fi

# unique outdir for this job; don't pre-create it
OUTDIR="${OUTDIR_BASE}_${SLURM_JOB_ID}"

# optional: node-local scratch
if [[ -d /scratch/$USER ]]; then
  export TMPDIR="/scratch/$USER/${SLURM_JOB_ID}"
  mkdir -p "$TMPDIR"
  echo "TMPDIR=$TMPDIR"
fi

# --- Roary (nt core alignment) ---
roary -e --mafft -p "${THREADS}" -i "${ID_PCT}" -cd "${CORE_PCT}" \
      -f "${OUTDIR}" -v "${GFFS[@]}"

# If Roary still added a suffix (paranoia), pick the newest roary_out*
if [[ ! -d "${OUTDIR}" ]]; then
  OUTDIR=$(ls -dt "${OUTDIR_BASE}"* 2>/dev/null | head -n1 || true)
fi

ALIGN="${OUTDIR}/core_gene_alignment.aln"
if [[ ! -s "${ALIGN}" ]]; then
  echo "ERROR: ${ALIGN} missing or empty"; ls -lh "${OUTDIR}" || true; exit 3
fi

# --- FastTree (nt, GTR+Gamma) ---
TREE="${OUTDIR}/core_gene_alignment.fasttree.nwk"
if command -v FastTreeMP >/dev/null 2>&1; then
  FastTreeMP -gtr -nt -gamma < "${ALIGN}" > "${TREE}"
else
  FastTree    -gtr -nt -gamma < "${ALIGN}" > "${TREE}"
fi

echo "Outputs:"
echo "  - ${OUTDIR}/gene_presence_absence.csv"
echo "  - ${ALIGN}"
echo "  - ${TREE}"
