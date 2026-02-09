
soft=$1
run=${PWD##*/}
nextflow  run  ${soft}/AD_analysis/main.nf \
        --batch ${run} \
        --input  data/ \
        --output ${PWD} \
        -entry sispa_workflow -profile standard -resume \
        -with-trace -with-report 