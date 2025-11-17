import warnings
# Silence all UserWarnings
warnings.filterwarnings('ignore', category=UserWarning)

import argparse
import sys
try:
    import pandas as pd
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
except ImportError as e:
    print(f"Error: {e}. Please make sure you have pandas and alphagenome installed. You can install them using pip: pip install pandas alphagenome-client", file=sys.stderr)
    sys.exit(1)

def score_variants(api_key, input_csv_path, species, output):
    """
    Scores variants from a CSV file using the AlphaGenome API.

    Args:
        api_key (str): Your AlphaGenome API key.
        input_csv_path (str): Path to the input CSV file.
        species (str): 'human' or 'mouse'.
        output (str): Path to the output CSV file.
    """
    dna_model = dna_client.create(api_key)
    try:
        df = pd.read_csv(input_csv_path)
    except FileNotFoundError:
        print(f"Error: The file '{input_csv_path}' was not found.")
        return

    # Map GRanges-style columns to VCF-style columns if necessary
    # GRanges -> VCF
    # seqnames -> CHROM
    # start -> POS
    column_mapping = {
        'seqnames': 'CHROM',
        'start': 'POS'
    }
    df = df.rename(columns=column_mapping)
    required_columns = ['CHROM', 'POS', 'REF', 'ALT']
    for col in required_columns:
        if col not in df.columns:
            print(f"Error: Input CSV is missing required column '{col}'.")
            return

    if 'variant_id' not in df.columns:
        df['variant_id'] = df.apply(
            lambda row: f"{row['CHROM']}:{row['POS']}:{row['REF']}>{row['ALT']}",
            axis=1
        )

    organism_map = {
        'human': dna_client.Organism.HOMO_SAPIENS,
        'mouse': dna_client.Organism.MUS_MUSCULUS,
    }
    organism = organism_map.get(species)
    if not organism:
        print(f"Error: Invalid species '{species}'. Please choose 'human' or 'mouse'.")
        return

    # Define the scorers for splicing and open chromatin
    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected_scorers = [
        all_scorers[key] for key in all_scorers
        if key.lower() in [
            'splice_sites', 'splice_site_usage', #'splice_junctions',
            'atac', 'dnase'
        ]
    ]
    
    # Remove unsupported scorers for the selected organism
    unsupported_scorers = [
        scorer for scorer in selected_scorers
        if organism.value not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
    ]
    if unsupported_scorers:
        for scorer in unsupported_scorers:
            selected_scorers.remove(scorer)
            
    # Process variants and collect scores
    results = []
    for _, row in df.iterrows():
        variant = genome.Variant(
            chromosome=str(row['CHROM']),
            position=int(row['POS']),
            reference_bases=row['REF'],
            alternate_bases=row['ALT'],
            name=row['variant_id'],
        )
        
        # The interval is the sequence length around the variant to be analyzed
        sequence_length = '16KB'  # @param ["16KB", "100KB", "500KB", "1MB"] { type:"string" }
        sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
            f'SEQUENCE_LENGTH_{sequence_length}'
        ]
        interval = variant.reference_interval.resize(sequence_length)
        variant_scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=selected_scorers,
            organism=organism,
        )
        results.append(variant_scores)

    df_scores = variant_scorers.tidy_scores(results)

    # Merge the scores with the original data
    # augmented_df = pd.merge(df, df_scores, on='variant_id', how='left')
    df_scores.to_csv(output, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Score variants using the AlphaGenome API.')
    parser.add_argument(
      '--api_key', type=str, required=True, help='Your AlphaGenome API key.')
    parser.add_argument(
      '--input', type=str, required=True, help='Path to the input CSV file with variant data.')
    parser.add_argument(
      '--species', type=str, required=True, 
      choices=['human', 'mouse'], help='Species for the analysis.')
    parser.add_argument(
      '--output', type=str, required=True, help='Output file CSV.')

    args = parser.parse_args()
    score_variants(args.api_key, args.input, args.species, args.output)
