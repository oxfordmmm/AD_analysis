# Run script from one level up. So command will be ./envs/create.sh

for e in envs/*.yaml
do 
    stem=$(basename $e .yaml)
    echo $stem
    mamba env create -f $e -p condaCache/$stem
done