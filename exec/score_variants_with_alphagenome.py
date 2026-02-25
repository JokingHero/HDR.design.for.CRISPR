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
    print(f"Error: {e}. Please make sure you have pandas and alphagenome installed.", file=sys.stderr)
    sys.exit(1)

def write_empty_output(output_path):
    """Writes an empty DataFrame with the required columns for R to parse."""
    empty_df = pd.DataFrame(columns=['variant_id', 'output_type', 'quantile_score'])
    empty_df.to_csv(output_path, index=False)

def score_variants(api_key, input_csv_path, species, output):
    """
    Scores variants from a CSV file using the AlphaGenome API.
    """
    try:
        dna_model = dna_client.create(api_key)
        try:
            df = pd.read_csv(input_csv_path)
        except FileNotFoundError:
            print(f"Error: The file '{input_csv_path}' was not found.", file=sys.stderr)
            write_empty_output(output)
            sys.exit(1)

        if df.empty:
            write_empty_output(output)
            return

        column_mapping = {
            'seqnames': 'CHROM',
            'start': 'POS'
        }
        df = df.rename(columns=column_mapping)
        required_columns = ['CHROM', 'POS', 'REF', 'ALT']
        for col in required_columns:
            if col not in df.columns:
                print(f"Error: Input CSV is missing required column '{col}'.", file=sys.stderr)
                write_empty_output(output)
                sys.exit(1)

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
            print(f"Error: Invalid species '{species}'. Please choose 'human' or 'mouse'.", file=sys.stderr)
            write_empty_output(output)
            sys.exit(1)

        all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
        selected_scorers = [
            all_scorers[key] for key in all_scorers
            if key.lower() in ['splice_sites', 'splice_site_usage', 'splice_junctions']
        ]
        
        unsupported_scorers = [
            scorer for scorer in selected_scorers
            if organism.value not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
        ]
        for scorer in unsupported_scorers:
            selected_scorers.remove(scorer)
                
        results = []
        for _, row in df.iterrows():
            variant = genome.Variant(
                chromosome=str(row['CHROM']),
                position=int(row['POS']),
                reference_bases=row['REF'],
                alternate_bases=row['ALT'],
                name=row['variant_id'],
            )
            
            sequence_length = '1MB'
            sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f'SEQUENCE_LENGTH_{sequence_length}']
            interval = variant.reference_interval.resize(sequence_length)
            
            try:
                variant_scores = dna_model.score_variant(
                    interval=interval,
                    variant=variant,
                    variant_scorers=selected_scorers,
                    organism=organism,
                )
                results.append(variant_scores)
            except Exception as e:
                print(f"Warning: Failed to score variant {row['variant_id']}: {e}", file=sys.stderr)

        if not results:
            print("Warning: No valid variants were scored.", file=sys.stderr)
            write_empty_output(output)
            return

        # Attempt to format the API results using Pandas
        try:
            df_scores = variant_scorers.tidy_scores(results)
        except KeyError as e:
            # THIS catches your specific error (if API returns empty payload)
            print(f"Warning: AlphaGenome API returned no usable scores (KeyError: {e}).", file=sys.stderr)
            write_empty_output(output)
            return

        df_scores.to_csv(output, index=False)

    except Exception as e:
        # Global catch-all to ensure R always receives an output file, even on total crashes.
        print(f"Critical Error in Python script: {e}", file=sys.stderr)
        write_empty_output(output)
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Score variants using the AlphaGenome API.')
    parser.add_argument('--api_key', type=str, required=True, help='Your AlphaGenome API key.')
    parser.add_argument('--input', type=str, required=True, help='Path to the input CSV file with variant data.')
    parser.add_argument('--species', type=str, required=True, choices=['human', 'mouse'], help='Species for the analysis.')
    parser.add_argument('--output', type=str, required=True, help='Output file CSV.')

    args = parser.parse_args()
    score_variants(args.api_key, args.input, args.species, args.output)
