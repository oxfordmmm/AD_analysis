# Agnostic Diagnostic

This is a pipeline for taking metagenomic reads amplified with Phi or Sispa. Currently focused on detecting 5 spikes: MS2, armored RNA, Zika, Murine Respirovirus, Orthoreovirus.

The process is based around detecting the primer (N9-GATGATAGTAGGGCTTCGTCAC) which is used in the reverse transcription step and by SISPA.

To replicate our validation study using this workflow, please see the [study replication folder](study_replication/README.md).

## SISPA
With SISPA pieces of DNA with primer bound to either end get amplified. So the expectation for these reads is that primer should be found at the ends of the ONT reads. These may not be found if the read was fragmented before going through the pore, or if the start/end aren't read so well.

### SISPA Process
1. Blastn: A tool called blastn is used to find all the places in the reads where the primer sequence appears.
2. Trimming: A script goes through every single read. If there is any primer near the start or end then it gets cutoff. If the read had no primer then it is left as it is.
3. mapping: A tool called minimap2 is used to align the reads to the references. For each read we only keep the best/primary alignment.
4. consensus: Now we have all the reads mapped we can form the consensus. Basically for each position in the reference you consider all the reads which map to that place and decide what base to call (potentially a SNP if enough evidence for being different to the reference). If the coverage at a position is too low (currently less than 5 reads), than we report it as ambiguous using the letter N.
5. Error counts: We can then look at the resulting consensus and count the number of N's (positions with not enough coverage) as well as the number of differences to the reference.
6. Coverage: The other thing you can do once you have all the reads aligned is to look at coverage/depth. For each position you can just count the number of reads which map to that position. I've got an R script to make a graph of this for each of the references.

## PHI
In PHI we form circles of DNA which should be up to 1kb long and contain the primer somewhere in the circle. PHI then loops round the circle amplifying it. This should generate long reads which should have the primer evenly spaced. However, these don't always go so cleanly and PHI may could and replicate one of the repeating reads already created and generally produce reads which are a bit messy.

### PHI Process
The process is similar to SISPA so only differences are described here.
1. Blastn: same as SISPA
2. Splitting: We go through every read. If it has no primer than leave it as is. Otherwise we split the read at every place where the primer is found to form several segments. Sometimes you find primer very close together (even back to back!), in these cases the resulting segments are very short so we just remove them (There is a min_contig_length parameter for the script). **Each of these segments will now be treated as a seperate read**. However their read id will have the added ending 'splitX' so we can keep track of what happens to all the different segments of a read.  
3. mapping: Again use minimap2. But this time we keep supplementary alignments (where different parts of the read map separately). This is really helpful. You find plenty of reads which we didn't detect the primer in yet still seem to have come from circles and so different parts of the read map to similar places in the reference. NOTE: we only keep alignments which don't overlap (too much).
4. Consensus: same as SISPA
5. Error counts: same as SISPA
6. Coverage: Same as SISPA. 
7. PHI stats: For each read we can try to get some idea of the circle it came from. For each read we can get some summary stats:
   1. Number of primers detecting in the read
   2. Number of segments this created. Note: because primers can be messy, this may be significantly lower than the number of primers.
   3. Number of alignments. Note: because some segments may give multiple alignments if they contain primers that were missed.
   4. Aligned bases. By summing the length of each alignment you can get a value for how many bases in total the read contributed.
   5. Circle length: By taking the median segment length you can get a sense of how large the original DNA circle must have been. Treat with caution.



## Outputs
Summary: In each output folder is a **summary** folder. This contains the most easy to digest outputs
1. error_rates.csv: This compares all barcodes for their percentage of bases which are N (coverage less than 5), and the percentage which are differences.
2. Nanoplots. This contains the plots of read length from the raw reads.
3. overall_summary.csv: This contains lots of general summary for each barcode. Most column names make sense, but here's some explanations.
   1. reads_all: This counts all the reads that make it through crumpit (and so get saved to rescomp).
   2. reads_mapped: this is the number of reads which mapped (potentially after splitting into pseudoreads) to any of the references. This is what mapped means in any of the columns
   3. have_primer: number of reads for which blast found the 22bp primer sequence (or large enough part of the 22bp)
   4. Q1, Q2, Q3 after the quartiles. Q95 is 95th percentile
   5. I'm labelling a read a 'pseudocircle' if it creates at least 3 pseudoreads when split by primer. The circle size is then the median length of its pseudoreads. If any of the pseudoreads are mapped than the whole circle is considered mapped.
   6. I'm labelling a read an align_circle if it has at least 3 non-overlapping alignments to the reference. These could come from different pseudoreads but not required. you can have a read with no primer found, but still has multiple alignments. The circle size is the median length of its alignments
   7. 'failed_trimming' is the sispa reads. The next three columns breakdown by reason for failing. Mostly this comes from having a primer somewhere in the middle of the read. But some had more than 10 primers which I also used as a filter.
 


---

1. Blastn. The **primer_blast** folder contains tsv which list every blast a primer was found.
2. Trimming/Splitting. 
   1. The **trimming** folder contains the trimmed reads, failed reads, and reads stats. There is also a txt report giving the number of reads in each category.
   2. The **splitting** folder contains info on each barcode. It has the split reads and a read stats csv. A report txt gives some more info on read counts.
3. mapping. The **bams** folder contains all the alignment files. minimap.bam is the output directly from minimap2. filt.bam is after the alignments have been filtered as described in the processes. 
Useful info about these alignments is put into the folder **alignment_info**. all_alignments.csv describes all the alignments minimap2 found. alignments.csv or non_overlapping.csv is after the alignments were filtered. Some python scripts process this to give a summary which you can find in **alignment_summary**
1. The **consensus** folder contains the consensus fastas. The **vcfs** contains the vcfs created by clair3.
2. folder **n_counts** contains files giving the N count and differences for each barcode.
3. Folder **depths** contains summaries of the average depth etc for each barcode. The **coverage_graphs** has the nice coverage graphs. The red points highlight places where coverage is less than 5. 
 


## Notes
1. Supplementary vs Secondary. For each read minimap will produce a **primary** alignment in which some part of the read maps to some part of the reference. A **secondary** alignment is an alternative place in the reference that the same part of the read could map to. A **supplementary** alignment is when a different part of the read maps somewhere in the reference. i.e.
   - bases 1-100 of the read map to 1500-1600 of zika. This is the primary alignment
   - bases 5-90 of the read map to 1400-1485 of zika. This is secondary
   - bases 120-210 of the read map to 1500-1590 of zika. This is supplementary.
   - bases 120-210 of the read map to 100-190 of zika. This is also supplementary.
2. Consensus and Clair3. Because ONT has a weird error profile we use a tool called clair3 in the consensus step. It's got some machine learning magic or something. It produces a vcf file which is just a list of the differences from the reference. These could be SNPs of INDELS.
3. Primer removal: The primer also has the random 9 bases. So whenever we split/trim based on a primer we actually also remove the 9 bases before/after it based on the direction of the primer. 
