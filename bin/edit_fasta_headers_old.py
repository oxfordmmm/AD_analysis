#!/usr/bin/env python3
import os

# Provide the path to your directory containing FASTA files
# If adding this to nextflow use the $refDir parameter and arg[1]
directory_path='/well/bag/users/vbr851/projects/agnostic_diagnostic_v2/references/refs'

def edit_fasta_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    edited_lines = []
    for line in lines:
        if line.startswith('>NC_'):
            parts = line.split(' ', 1)
            if len(parts) > 1:
                first_word_removed = parts[1].strip()
                joined_words = ', '.join(['_'.join(word.strip().split(' ')) for word in first_word_removed.split(',')])
                joined_words = joined_words.replace(',', '').replace('(', '').replace(')', '').replace('/','_')
                edited_header = '>' + joined_words
                edited_lines.append(edited_header + '\n')
        else:
            edited_lines.append(line)

    with open(file_path, 'w') as file:
        file.writelines(edited_lines)

for filename in os.listdir(directory_path):
    if filename.endswith(".fasta"):
        edit_fasta_file(os.path.join(directory_path, filename))