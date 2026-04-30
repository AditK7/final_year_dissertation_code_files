#!/bin/bash
# ============================================
# SLURM Array Job - AutoDock Vina HTVS
# Docks all prepared ligands in parallel
# ============================================

#SBATCH --job-name=vina_htvs
#SBATCH --partition=defq
#SBATCH --array=1-1000
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2g
#SBATCH --time=1:00:00
#SBATCH --output=slurm_logs/vina_%A_%a.out
#SBATCH --error=slurm_logs/vina_%A_%a.err

# ============================================
# SETUP
# ============================================
module load autodockvina-uon/binary/1.2.5

WORK_DIR="${SLURM_SUBMIT_DIR}"
LIGAND_DIR="${WORK_DIR}/ligands_pdbqt"
RESULTS_DIR="${WORK_DIR}/htvs_results"
RECEPTOR="${WORK_DIR}/3Q07_rec.pdbqt"

# Get ligand number from array task ID
LIGAND_NUM=${SLURM_ARRAY_TASK_ID}
LIGAND_IN="${LIGAND_DIR}/ligand_${LIGAND_NUM}.pdbqt"
LIGAND_OUT="${RESULTS_DIR}/ligand_${LIGAND_NUM}_docked.pdbqt"
LOG_FILE="${RESULTS_DIR}/ligand_${LIGAND_NUM}_log.txt"

mkdir -p "$RESULTS_DIR"

# Skip if ligand doesn't exist (some may have failed prep)
if [ ! -f "$LIGAND_IN" ]; then
    echo "Ligand $LIGAND_NUM not found, skipping."
    exit 0
fi

# Skip if already docked (in case of resubmission)
if [ -f "$LIGAND_OUT" ]; then
    echo "Ligand $LIGAND_NUM already docked, skipping."
    exit 0
fi

# ============================================
# DOCK
# ============================================
echo "Docking ligand $LIGAND_NUM on $(hostname) at $(date)"

vina_1.2.5_linux_x86_64 \
    --receptor "$RECEPTOR" \
    --ligand "$LIGAND_IN" \
    --center_x -7.937667 \
    --center_y -16.306053 \
    --center_z 11.782777 \
    --size_x 17.755 \
    --size_y 13.701 \
    --size_z 10.238 \
    --exhaustiveness 64 \
    --num_modes 20 \
    --energy_range 4 \
    --out "$LIGAND_OUT" \
    > "$LOG_FILE" 2>&1

# Extract and print best score
BEST_SCORE=$(grep "^   1 " "$LOG_FILE" | awk '{print $2}')
echo "Ligand $LIGAND_NUM: Best score = $BEST_SCORE kcal/mol"
echo "Finished at $(date)"
