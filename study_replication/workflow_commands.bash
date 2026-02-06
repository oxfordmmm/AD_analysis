
soft=$1
run=${PWD##*/}
nextflow  run  ${soft}/AD_analysis/main.nf \
        --batch ${run} \
        --input  data/ \
        --output /mnt/data/analysis/nick/agnostic_diagnostic/AD_winter_study/${run} \
        -entry sispa_workflow -profile standard -resume \
        -with-trace -with-report 