runs="AD_winter_study_201224
AD_winter_study_220125
AD_winter_study_070325
AD_winter_study_170325
AD_winter_study_100425
AD_winter_study_090525
AD_winter_study_160725
AD_winter_study_220725
AD_winter_study_240725
AD_winter_study_300725
AD_winter_study_010825
AD_winter_study_130825_rpt050825"

soft=$1

for run in ${runs}
do
        mkdir -p ${run}
        cp workflow_commands.bash ${run}/
        cd ${run}
        bash workflow_commands.bash ${soft}
        cd ..
done
