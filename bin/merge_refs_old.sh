#!/bin/bash

refDir=/well/bag/users/vbr851/projects/agnostic_diagnostic_v2/references
echo "Creating metareference in $refDir/"
cat $refDir/refs/*.fasta > $refDir/meta_ref.fasta