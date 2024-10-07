#!/bin/bash

# Input file path
input=$1

echo "$input" | sed -e 's/,/ /g' -e 's/\(([^ ]\)/_/g' -e 's/[(|)]//g' -e 's/\//_/g' -e 's/ /_/g' -e 's/,/_/g' -e 's/_\+/_/g' 

