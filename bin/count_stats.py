import sys
import pandas as pd

def calculate_statistics(csv_file, output_file):
    try:
        # Read the CSV file into a pandas DataFrame
        df = pd.read_csv(csv_file)

        # Filter rows where 'has_primer' is true
        has_primer_true_df = df[df['has_primer'] == True]

        # Count the number of True and False values in the "has_primer" column
        true_count = has_primer_true_df['has_primer'].sum()
        false_count = len(df) - true_count

        # Count total number of reads
        total_read_sum = len(df)

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

        # Print the statistics to the console
        print(f"Total number of reads: {total_read_sum}")
        print(f"Number of reads with primers found: {true_count}")
        print(f"Number of non-primer reads: {false_count}")
        print(f"Total number of pseudoreads after splitting: {pseudoread_sum}")
        print(f"Mean number of pseudoreads for reads that were split: {pseudoread_mean:.2f}")
        print(f"Median number of pseudoreads for reads that were split: {pseudoread_median}")
        print(f"Mean of primer_count for reads that were split: {primer_count_mean:.2f}")
        print(f"Median of primer_count for reads that were split: {primer_count_median}")

        # Save the statistics to the output text file
        with open(output_file, 'w') as f:
            f.write(f"Total number of reads: {total_read_sum}\n")
            f.write(f"Number of reads with primers found: {true_count}\n")
            f.write(f"Number of non-primer reads: {false_count}\n")
            f.write(f"Total number of pseudoreads after splitting: {pseudoread_sum}\n")
            f.write(f"Mean number of pseudoreads for reads that were split: {pseudoread_mean:.2f}\n")
            f.write(f"Median number of pseudoreads for reads that were split: {pseudoread_median}\n")
            f.write(f"Mean of primer_count for reads that were split: {primer_count_mean:.2f}\n")
            f.write(f"Median of primer_count for reads that were split: {primer_count_median}\n")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_statistics.py <csv_file> <output_file>")
        sys.exit(1)

    csv_file = sys.argv[1]
    output_file = sys.argv[2]
    calculate_statistics(csv_file, output_file)



