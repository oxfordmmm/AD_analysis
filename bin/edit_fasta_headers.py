#!/usr/bin/env python3
import sys

# Read in fasta file, edit headers, and write to a new file
def edit_fasta_file(input_file_path, output_file_path):
    with open(input_file_path, 'r') as file:
        lines = file.readlines()

    edited_lines = []
    for line in lines:
        if line.startswith('>NC_'):
            parts = line.split(' ', 1)
            if len(parts) > 1:
                first_word_removed = parts[1].strip()
                joined_words = ', '.join(['_'.join(word.strip().split(' ')) for word in first_word_removed.split(',')])
                joined_words = joined_words.replace(',', '').replace('(', '').replace(')', '').replace('/', '_')
                edited_header = '>' + joined_words
                edited_lines.append(edited_header + '\n')
        else:
            edited_lines.append(line)

    with open(output_file_path, 'w') as file:
        file.writelines(edited_lines)

if __name__ == "__main__":
    # Check if both input and output file paths are provided as command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file_path> <output_file_path>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    edit_fasta_file(input_file_path, output_file_path)
