import sys
import pandas as pd

def calculate_statistics(csv_file):
    try:
        # Read the CSV file into a pandas DataFrame
        df = pd.read_csv(csv_file)

        # Filter rows where 'has_primer' is true
        has_primer_true_df = df[df['has_primer'] == True]

        # Count the number of True and False values in the "has_primer" column
        true_count = has_primer_true_df['has_primer'].sum()
        false_count = len(df) - true_count

        # Sum the values in the "pseudoread_count" column
        pseudoread_sum = has_primer_true_df['pseudoread_count'].sum()

        # Calculate the mean of the "pseudoread_count" column for rows with 'has_primer' true
        pseudoread_mean = has_primer_true_df['pseudoread_count'].mean()

        # Calculate the median of the "pseudoread_count" column for rows with 'has_primer' true
        pseudoread_median = has_primer_true_df['pseudoread_count'].median()

        # Calculate the mean of the "primer_count" column for rows with 'has_primer' true
        primer_count_mean = has_primer_true_df['primer_count'].mean()

        # Calculate the median of the "primer_count" column for rows with 'has_primer' true
        primer_count_median = has_primer_true_df['primer_count'].median()

        print(f"Number of True values: {true_count}")
        print(f"Number of False values: {false_count}")
        print(f"Total number of pseudoreads: {pseudoread_sum}")
        print(f"Mean of pseudoread_count (where 'has_primer' is true): {pseudoread_mean:.2f}")
        print(f"Median of pseudoread_count (where 'has_primer' is true): {pseudoread_median}")
        print(f"Mean of primer_count (where 'has_primer' is true): {primer_count_mean:.2f}")
        print(f"Median of primer_count (where 'has_primer' is true): {primer_count_median}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python calculate_statistics.py <csv_file>")
        sys.exit(1)

    csv_file = sys.argv[1]
    calculate_statistics(csv_file)

