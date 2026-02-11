runs="AD_winter_study_201224
AD_winter_study_220125
AD_winter_study_070325
AD_winter_study_170325
AD_winter_study_100425
AD_winter_study_160725
AD_winter_study_220725
AD_winter_study_240725
AD_winter_study_300725
AD_winter_study_010825
AD_winter_study_130825_rpt050825"

soft=$1

baseDir=${PWD}

for run in ${runs}
do
        #mkdir -p ${run}
        cp workflow_commands.bash workflow_runs/${run}/
        cd workflow_runs/${run}
        bash workflow_commands.bash ${soft}
        cd ${baseDir}
done
