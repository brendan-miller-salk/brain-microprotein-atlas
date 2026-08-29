import os
import pandas as pd

# Define the root directory
root_dir = "/scratch1/brendajm/tmt_rosmap"

# Initialize a list to hold the combined data from all subdirectories
combined_data = []

# Loop through each subdirectory with "b" prefix
for i in range(1, 51):  # Assuming directories are b1 to b50
    sub_dir = os.path.join(root_dir, f"b{i}/shortstop_proteogenomics_appended_results_cpm05/DDA")
    protein_file = os.path.join(sub_dir, "protein.tsv")

    # Check if the file exists
    if os.path.exists(protein_file):
        # Read the protein.tsv file
        df = pd.read_csv(protein_file, sep="\t")

        # Extract relevant columns
        if all(col in df.columns for col in ['Protein', 'Unique Spectral Count', 'Razor Spectral Count', 'Indistinguishable Proteins']):
            extracted_data = df[['Protein', 'Unique Spectral Count', 'Razor Spectral Count', 'Indistinguishable Proteins']]

            # Append the data to the combined list
            combined_data.append(extracted_data)
        else:
            print(f"Missing required columns in file: {protein_file}")
    else:
        print(f"File not found: {protein_file}")

# Combine all the data into a single DataFrame
if combined_data:
    combined_df = pd.concat(combined_data, ignore_index=True)

    # Replace NaN or empty Indistinguishable Proteins with "BLANK"
    combined_df['Indistinguishable Proteins'] = combined_df['Indistinguishable Proteins'].fillna("BLANK")
    combined_df.loc[combined_df['Indistinguishable Proteins'] == "", 'Indistinguishable Proteins'] = "BLANK"

    # Group by Protein and aggregate the data
    aggregated_df = combined_df.groupby(
        'Protein', as_index=False
    ).agg({
        'Unique Spectral Count': 'sum',
        'Razor Spectral Count': 'sum',
        'Indistinguishable Proteins': lambda x: ";".join(set(x))
    })

    # Save the aggregated data to a CSV file
    aggregated_csv = os.path.join(root_dir, "shortstop_proteogenomics_appended_proteinID_uniqueness.csv")
    aggregated_df.to_csv(aggregated_csv, index=False)

    # Save the Protein IDs to a text file
    protein_ids = os.path.join(root_dir, "shortstop_proteogenomics_appended_proteinID_uniqueness.txt")
    aggregated_df['Protein'].to_csv(protein_ids, index=False)

    print(f"Aggregated data saved to {aggregated_csv}")
else:
    print("No data found across the directories.")
