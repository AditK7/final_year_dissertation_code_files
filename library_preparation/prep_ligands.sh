#!/bin/bash
# ============================================
# Ligand Preparation Script
# Splits SDF and converts all ligands to PDBQT
# Run this BEFORE submitting the SLURM docking job
# ============================================

SDF_FILE="chembl_1000_compounds.sdf"
LIGAND_DIR="ligands_pdbqt"

mkdir -p "$LIGAND_DIR"

# Count compounds
TOTAL=$(grep -c '$$$$' "$SDF_FILE")
echo "Found $TOTAL compounds in $SDF_FILE"
echo "Splitting and converting to PDBQT..."
echo ""

SUCCESS=0
FAIL=0

for i in $(seq 1 $TOTAL); do
    SDF_OUT="$LIGAND_DIR/ligand_${i}.sdf"
    MOL2_OUT="$LIGAND_DIR/ligand_${i}.mol2"
    PDBQT_OUT="$LIGAND_DIR/ligand_${i}.pdbqt"

    # Split
    obabel "$SDF_FILE" -O "$SDF_OUT" -f $i -l $i 2>/dev/null

    # SDF -> MOL2 with 3D coords and hydrogens
    obabel "$SDF_OUT" -O "$MOL2_OUT" --gen3d -h 2>/dev/null
    if [ $? -ne 0 ] || [ ! -s "$MOL2_OUT" ]; then
        echo "  FAILED prep: ligand_${i} (obabel)"
        FAIL=$((FAIL + 1))
        continue
    fi

    # MOL2 -> PDBQT with Meeko
    mk_prepare_ligand.py -i "$MOL2_OUT" -o "$PDBQT_OUT" 2>/dev/null
    if [ $? -ne 0 ] || [ ! -s "$PDBQT_OUT" ]; then
        echo "  FAILED prep: ligand_${i} (meeko)"
        FAIL=$((FAIL + 1))
        continue
    fi

    SUCCESS=$((SUCCESS + 1))

    # Progress update every 50
    if [ $((i % 50)) -eq 0 ]; then
        echo "  Processed $i / $TOTAL ($SUCCESS success, $FAIL failed)"
    fi
done

echo ""
echo "=========================================="
echo "LIGAND PREPARATION COMPLETE"
echo "=========================================="
echo "Total:      $TOTAL"
echo "Success:    $SUCCESS"
echo "Failed:     $FAIL"
echo "PDBQT dir:  $LIGAND_DIR/"

# Write list of successfully prepared ligands
echo ""
echo "Writing ligand list..."
ls "$LIGAND_DIR"/*.pdbqt | grep -v "_docked" | wc -l > "$LIGAND_DIR/total_count.txt"
echo "Ready for SLURM submission."
