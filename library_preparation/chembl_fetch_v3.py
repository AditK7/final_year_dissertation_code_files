#!/usr/bin/env python3
"""
ChEMBL Compound Fetcher v3 - HPC Edition
Fetches 1000 compounds in MW range 467-567 for HTVS
"""

import requests
import json
import time
import re
import sys

# ============================================
# SETTINGS
# ============================================
MW_MIN = 467
MW_MAX = 567
MAX_COMPOUNDS = 1000
ALLOWED_ATOMS = {'H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'}
OUTPUT_SDF = "chembl_1000_compounds.sdf"
OUTPUT_LOG = "chembl_1000_compound_ids.csv"

API_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"

print("=" * 65)
print("ChEMBL Compound Fetcher v3 - HPC Edition")
print("=" * 65)
print(f"\nCriteria:")
print(f"  Molecular weight: {MW_MIN} - {MW_MAX} Da")
print(f"  Molecule type:    Small molecule")
print(f"  Allowed atoms:    {', '.join(sorted(ALLOWED_ATOMS))}")
print(f"  Target output:    {MAX_COMPOUNDS} compounds\n")

# Step 1: Get total count
print("[1/5] Querying total matching compounds in ChEMBL...")
count_params = {
    "molecule_properties__full_mwt__gte": MW_MIN,
    "molecule_properties__full_mwt__lte": MW_MAX,
    "molecule_type": "Small molecule",
    "limit": 1,
}
try:
    resp = requests.get(API_URL, params=count_params, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    total_in_chembl = data.get("page_meta", {}).get("total_count", "unknown")
    print(f"       Total in ChEMBL: {total_in_chembl}")
except Exception as e:
    print(f"       Could not get total count: {e}")
    total_in_chembl = "unknown"

# Step 2: Fetch compounds
print(f"\n[2/5] Fetching compounds...")
params = {
    "molecule_properties__full_mwt__gte": MW_MIN,
    "molecule_properties__full_mwt__lte": MW_MAX,
    "molecule_type": "Small molecule",
    "limit": 1000,
    "offset": 0,
}

all_compounds = []
batch = 0
fetch_target = MAX_COMPOUNDS * 3

while len(all_compounds) < fetch_target:
    batch += 1
    print(f"       Batch {batch}: offset {params['offset']} "
          f"({len(all_compounds)} so far)...")
    try:
        resp = requests.get(API_URL, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        print(f"       Error: {e}")
        break

    molecules = data.get("molecules", [])
    if not molecules:
        break

    for mol in molecules:
        chembl_id = mol.get("molecule_chembl_id", "")
        structs = mol.get("molecule_structures")
        props = mol.get("molecule_properties")
        if not structs or not props:
            continue
        smiles = structs.get("canonical_smiles", "")
        molfile = structs.get("molfile")
        mw = props.get("full_mwt")
        formula = props.get("molecular_formula", "")
        pref_name = mol.get("pref_name", "N/A") or "N/A"
        if not smiles or not molfile:
            continue
        all_compounds.append({
            "chembl_id": chembl_id,
            "smiles": smiles,
            "molfile": molfile,
            "mw": mw,
            "formula": formula,
            "pref_name": pref_name,
        })

    if not data.get("page_meta", {}).get("next"):
        break
    params["offset"] += 1000
    time.sleep(0.5)

print(f"       Fetched {len(all_compounds)} compounds with structures.")

# Step 3: Filter exotic atoms
print(f"\n[3/5] Filtering exotic atoms...")

def has_only_allowed_atoms(formula):
    elements = re.findall(r'[A-Z][a-z]?', formula)
    return all(el in ALLOWED_ATOMS for el in elements)

clean_compounds = [c for c in all_compounds if has_only_allowed_atoms(c["formula"])]
print(f"       Passed: {len(clean_compounds)}")
print(f"       Removed: {len(all_compounds) - len(clean_compounds)}")

# Step 4: Select
selected = clean_compounds[:MAX_COMPOUNDS]
print(f"\n[4/5] Selected {len(selected)} compounds.")

# Step 5: Write outputs
print(f"\n[5/5] Writing output files...")

with open(OUTPUT_SDF, "w") as f:
    for i, comp in enumerate(selected, 1):
        molfile = comp["molfile"]
        f.write(molfile)
        if not molfile.strip().endswith("M  END"):
            f.write("M  END\n")
        f.write(f">  <chembl_id>\n{comp['chembl_id']}\n\n")
        f.write(f">  <ligand_number>\n{i}\n\n")
        f.write(f">  <smiles>\n{comp['smiles']}\n\n")
        f.write(f">  <mw>\n{comp['mw']}\n\n")
        f.write(f">  <formula>\n{comp['formula']}\n\n")
        f.write(f">  <pref_name>\n{comp['pref_name']}\n\n")
        f.write("$$$$\n")

print(f"       SDF saved: {OUTPUT_SDF}")

with open(OUTPUT_LOG, "w") as f:
    f.write("ligand_number,chembl_id,pref_name,molecular_weight,formula,smiles\n")
    for i, comp in enumerate(selected, 1):
        smiles_clean = comp['smiles'].replace(',', ';')
        f.write(f"{i},{comp['chembl_id']},{comp['pref_name']},"
                f"{comp['mw']},{comp['formula']},{smiles_clean}\n")

print(f"       CSV saved: {OUTPUT_LOG}")

print(f"\n{'=' * 65}")
print(f"SUMMARY")
print(f"{'=' * 65}")
print(f"  Total in ChEMBL (MW {MW_MIN}-{MW_MAX}): {total_in_chembl}")
print(f"  Fetched with structures:   {len(all_compounds)}")
print(f"  Passed atom filter:        {len(clean_compounds)}")
print(f"  Selected for docking:      {len(selected)}")
print(f"  Output SDF:                {OUTPUT_SDF}")
print(f"  Output CSV:                {OUTPUT_LOG}")
