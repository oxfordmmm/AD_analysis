# Agnostic Diagnostic v2

Currently, there are two main nextflow workflows:
1. pre_main.nf - This downloads all of the references based on a taxid list and creates a meta-reference
2. main.nf

- references/refs already contains the references (.fasta and .gff) for the control viruses we spiked in
- references/taxid_list.txt contains the list of taxids for the references we want to download (not including the controls)
- the first workflow will download these extra refs and create a metaref.fasta file containing all of the references we want to map to
- do we actually need the gff files?
- in the future, add the taxon ids for the controls to the taxid_list.txt OR if these fastq files are not in the refseq database and metadata then create a text file with a list of virus names so that we can automate the plotting scripts later - to group any fastq headers as fastq file names for the mapping bar plots