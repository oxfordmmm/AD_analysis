#!/bin/bash

fastas=$1

echo "Creating meta-reference"

#echo "Creating metareference in $refDir/"
#cat $refDir/refs/*.fasta > $refDir/meta_ref.fasta

cat $fastas > meta_ref.fasta