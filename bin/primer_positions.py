#!/usr/bin/env python3
# Define a dictionary to store sseqid as keys and a list of (sstart, send) positions as values
primer_positions = {}

# Read the tabular data line by line (excluding the header line)
with open('1_blast.tsv', 'r') as file:
    next(file)  # Skip the header line
    for line in file:
        fields = line.strip().split('\t')
        sseqid = fields[1]
        sstart = int(fields[8])
        send = int(fields[9])
        
        # Check if the sseqid is already in the dictionary
        if sseqid in primer_positions:
            primer_positions[sseqid].append((sstart, send))
        else:
            primer_positions[sseqid] = [(sstart, send)]

# Now, primer_positions is a dictionary where keys are sseqid and values are lists of (sstart, send) positions

# Create a new TSV file to write the primer positions
with open('primer_positions.tsv', 'w') as output_file:
    output_file.write("sseqid\tprimer_positions\n")
    
    for sseqid, positions in primer_positions.items():
        primer_str = ", ".join([f"({start}, {end})" for start, end in positions])
        output_file.write(f"{sseqid}\t[{primer_str}]\n")

