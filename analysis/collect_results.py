#!/usr/bin/env python3
"""
Collect and rank HTVS docking results.
Run after all SLURM array jobs have completed.
"""

import os
import re
import csv

RESULTS_DIR = "htvs_results"
COMPOUND_LOG = "chembl_1000_compound_ids.csv"
OUTPUT_FILE = "htvs_ranked_results.csv"

# Load ChEMBL ID mapping
chembl_map = {}
if os.path.exists(COMPOUND_LOG):
    with open(COMPOUND_LOG, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            num = int(row['ligand_number'])
            chembl_map[num] = {
                'chembl_id': row['chembl_id'],
                'pref_name': row['pref_name'],
                'mw': row['molecular_weight'],
                'formula': row['formula'],
            }

# Collect results
results = []
for f in os.listdir(RESULTS_DIR):
    if f.endswith("_log.txt"):
        match = re.match(r'ligand_(\d+)_log\.txt', f)
        if not match:
            continue
        num = int(match.group(1))
        filepath = os.path.join(RESULTS_DIR, f)

        with open(filepath, 'r') as fh:
            content = fh.read()

        # Extract best score (mode 1)
        score_match = re.search(r'^\s+1\s+([-\d.]+)', content, re.MULTILINE)
        if score_match:
            score = float(score_match.group(1))
            info = chembl_map.get(num, {})
            results.append({
                'ligand_num': num,
                'score': score,
                'chembl_id': info.get('chembl_id', 'N/A'),
                'pref_name': info.get('pref_name', 'N/A'),
                'mw': info.get('mw', 'N/A'),
                'formula': info.get('formula', 'N/A'),
            })

# Sort by score (most negative = best)
results.sort(key=lambda x: x['score'])

# Print ranked results
print("=" * 90)
print("HTVS RESULTS - RANKED BY BINDING AFFINITY")
print("=" * 90)
print(f"\nTotal successfully docked: {len(results)}")
print(f"\n{'Rank':<6} {'Ligand':<10} {'Score':<12} {'ChEMBL ID':<16} {'MW':<10} {'Name'}")
print(f"{'-'*6} {'-'*9} {'-'*11} {'-'*15} {'-'*9} {'-'*20}")

for rank, r in enumerate(results, 1):
    name = r['pref_name'][:25]
    print(f"{rank:<6} L_{r['ligand_num']:<7} {r['score']:<12.3f} "
          f"{r['chembl_id']:<16} {r['mw']:<10} {name}")

    # Show top 20 in detail, then just top line
    if rank == 20:
        print(f"\n... showing top 20 of {len(results)} total ...")
        break

# Save full results to CSV
with open(OUTPUT_FILE, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['rank', 'ligand_number', 'vina_score_kcal_mol',
                      'chembl_id', 'pref_name', 'molecular_weight', 'formula'])
    for rank, r in enumerate(results, 1):
        writer.writerow([rank, r['ligand_num'], r['score'],
                         r['chembl_id'], r['pref_name'], r['mw'], r['formula']])

print(f"\nFull ranked results saved to: {OUTPUT_FILE}")

# Summary stats
scores = [r['score'] for r in results]
print(f"\n{'=' * 90}")
print(f"SUMMARY STATISTICS")
print(f"{'=' * 90}")
print(f"  Total docked:    {len(scores)}")
print(f"  Best score:      {min(scores):.3f} kcal/mol")
print(f"  Worst score:     {max(scores):.3f} kcal/mol")
print(f"  Mean score:      {sum(scores)/len(scores):.3f} kcal/mol")
print(f"  Scores < -8.0:   {sum(1 for s in scores if s < -8.0)}")
print(f"  Scores < -7.0:   {sum(1 for s in scores if s < -7.0)}")
print(f"  Scores < -6.0:   {sum(1 for s in scores if s < -6.0)}")
