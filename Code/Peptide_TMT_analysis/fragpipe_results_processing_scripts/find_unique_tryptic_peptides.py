#!/usr/bin/env python3
import os
import re
import sys
import argparse
import pandas as pd
import concurrent.futures

def process_peptide_file(peptide_file, gene_ids):
    """
    Process a single peptide file and return a list of dicts with gene_id and peptide.
    Uses a combined regex for efficiency.
    """
    results = []
    try:
        df = pd.read_csv(peptide_file, sep="\t", dtype=str)
    except Exception as e:
        print(f"Error reading {peptide_file}: {e}")
        return results

    # Check for required columns
    if "Protein" not in df.columns or "Peptide" not in df.columns:
        print(f"Missing required columns in {peptide_file}")
        return results

    # Exclude rows that contain "|" in any column
    mask_no_pipe = ~df.apply(lambda row: row.astype(str).str.contains("\|").any(), axis=1)
    df = df[mask_no_pipe]

    if df.empty:
        return results

    # Build a regex pattern that matches any gene id.
    # Using re.escape ensures that special characters in gene ids are handled properly.
    pattern = '|'.join(re.escape(g) for g in gene_ids)

    # Filter rows where Protein contains any of the gene ids (vectorized)
    matched_df = df[df["Protein"].str.contains(pattern, na=False)]
    if matched_df.empty:
        return results

    # For each row, use re.findall to determine which gene ids match.
    for _, row in matched_df.iterrows():
        protein_str = row["Protein"]
        peptide = row["Peptide"]
        peptide_start = row["Start"]
        peptide_end = row["End"]
        matching_genes = re.findall(pattern, protein_str)
        # Remove duplicate gene ids (if multiple matches)
        matching_genes = list(set(matching_genes))
        for gene in matching_genes:
            results.append({"gene_id": gene, "tryptic_peptide": peptide, "tryptic_peptide_start": peptide_start, "tryptic_peptide_end": peptide_end})
    return results

def find_unique_microprotein_peptides_bulk(gene_ids_file, main_dir, output_csv, num_workers=4):
    """
    Process directories b1 to b50 under main_dir, search peptide.tsv files,
    and save gene_id-peptide matches into output_csv using parallel processing.
    """
    # Read gene IDs (one per line)
    try:
        with open(gene_ids_file, 'r') as f:
            gene_ids = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Error reading gene IDs from {gene_ids_file}: {e}")
        return pd.DataFrame()

    print(f"Gene IDs: {gene_ids}")
    all_results = []
    peptide_files = []

    # Build list of peptide.tsv file paths from subdirectories b1 to b50
    for i in range(1, 51):
        sub_dir = os.path.join(main_dir, f"b{i}", "shortstop_proteogenomics_appended_results_cpm05", "DDA")
        peptide_file = os.path.join(sub_dir, "peptide.tsv")
        if os.path.exists(peptide_file):
            peptide_files.append(peptide_file)
        else:
            print(f"File not found: {peptide_file}")

    for i in range(1, 14):
        sub_dir = os.path.join(main_dir, f"round2/b{i}", "shortstop_proteogenomics_appended_results_cpm05", "DDA")
        peptide_file = os.path.join(sub_dir, "peptide.tsv")
        if os.path.exists(peptide_file):
            peptide_files.append(peptide_file)
        else:
            print(f"File not found: {peptide_file}")

    # Process files in parallel using ProcessPoolExecutor
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit a future for each peptide file
        future_to_file = {executor.submit(process_peptide_file, pf, gene_ids): pf for pf in peptide_files}
        for future in concurrent.futures.as_completed(future_to_file):
            pf = future_to_file[future]
            try:
                results = future.result()
                print(f"Found {len(results)} matches in {pf}")
                all_results.extend(results)
            except Exception as e:
                print(f"Error processing {pf}: {e}")

    if all_results:
        out_df = pd.DataFrame(all_results)
        out_df.drop_duplicates(inplace=True)
        out_df.to_csv(output_csv, index=False)
        print(f"Finished! Results saved in {output_csv}")
        return out_df
    else:
        print("No matches found.")
        return pd.DataFrame()

def main():
    parser = argparse.ArgumentParser(description="Find unique microprotein peptides with parallel processing.")
    parser.add_argument("gene_ids_file", help="File with gene IDs (one per line)")
    parser.add_argument("main_dir", help="Main directory containing subdirectories b1 to b50")
    parser.add_argument("output_csv", help="Output CSV file to save results")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers (default: 4)")
    args = parser.parse_args()

    df_out = find_unique_microprotein_peptides_bulk(args.gene_ids_file, args.main_dir, args.output_csv, num_workers=args.workers)
    print(df_out.head())

if __name__ == '__main__':
    main()
