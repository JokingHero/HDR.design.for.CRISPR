import argparse
import sys
import pandas as pd
import json
from alphagenome.models import dna_client

def fetch_splicing_metadata(api_key, output_csv):
    """
    Fetches metadata for all Human and Mouse tracks, filters for SPLICING,
    and saves a CSV including the 'species' column.
    """
    try:
        dna_model = dna_client.create(api_key)
    except Exception as e:
        print(f"Error creating client: {e}", file=sys.stderr)
        sys.exit(1)

    all_metadata = []

    # Define species mapping
    organisms = [
        ('human', dna_client.Organism.HOMO_SAPIENS),
        ('mouse', dna_client.Organism.MUS_MUSCULUS)
    ]

    for species_name, organism_enum in organisms:
        try:
            # Retrieve metadata dataframe
            meta_df = dna_model.output_metadata(organism=organism_enum).concatenate()
            meta_df['species'] = species_name
            all_metadata.append(meta_df)
            
        except Exception as e:
            print(f"Error fetching {species_name}: {e}", file=sys.stderr)

    if not all_metadata:
        print("No metadata retrieved.", file=sys.stderr)
        sys.exit(1)

    # Combine Human and Mouse into one DataFrame
    full_df = pd.concat(all_metadata, ignore_index=True)

    # --- FILTERING ---
    
    # 1. Filter for Splicing Only (Usage, Sites, Donor, Acceptor)
    target_types = ['OutputType.SPLICE_SITES', 'OutputType.SPLICE_SITE_USAGE']
    is_splicing = full_df['output_type'].astype(str).isin(target_types)
    splice_df = full_df[is_splicing].copy()

    # 2. Select Columns (Keeping 'species' as the first column)
    cols_to_keep = [
        'species', 
        'biosample_name'
    ]
    
    # Safety check to ensure columns exist
    available_cols = [c for c in cols_to_keep if c in splice_df.columns]
    final_df = splice_df[available_cols]
    final_df = final_df.sort_values(by=['species', 'biosample_name'])
    final_df = final_df.dropna(subset=['biosample_name'])
    final_df = final_df.drop_duplicates(subset=['species', 'biosample_name'])
    grouped_data = final_df.groupby('species')['biosample_name'].apply(list)
    output_dict = grouped_data.to_dict()
    with open(output, 'w') as f:
        # indent=2 makes it human-readable (pretty printed)
        json.dump(output_dict, f, indent=2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract AlphaGenome metadata for R package.')
    parser.add_argument('--api_key', type=str, required=True, help='Your AlphaGenome API key.')
    parser.add_argument('--output', type=str, required=True, help='Path to output CSV.')

    args = parser.parse_args()
    fetch_splicing_metadata(args.api_key, args.output)
