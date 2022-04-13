#!/bin/zsh
PLEIADESROOT=/nobackup/eanders/
TRUEPLEIADESROOT=$PLEIADESROOT/Empirical-Magnetic-Thermohaline/mesa_models/post/
FOLDERBASE=work_thermohaline

LOCALBASE=./
for tag in Brown_Garaud_Stellmach_13_coeff1.0e+00 
do
    for Zv in 0.00068468 0.00108379 0.00171429 0.00270849 0.00427157 0.00671765 0.01051793 0.01635613 0.02517204 0.03814436
    do
        for Mv in 0.9 1.1 1.3 1.5 1.7
        do
            echo $Zv $Mv
            FOLDEREND=alpha1.6000_Z$Zv\_M$Mv
            LOCALFOLDER=$FOLDERBASE\_$tag\_$FOLDEREND
            LOCALDIR=$LOCALBASE/$LOCALFOLDER
            PLEIADESDIR=$TRUEPLEIADESROOT/$LOCALFOLDER
            mkdir $LOCALDIR/
            files={figs_R0.mp4,r_vs_time.h5}
            scp pfe:$PLEIADESDIR/$files ./$LOCALDIR/
#            scp pfe:$PLEIADESDIR/fluid_parameter_plots.mp4 ./$LOCALDIR/
            figs={full_log_g_vs_model.png,full_r_vs_model.png,full_r_vs_log_g.png,full_R0_vs_model.png}
            scp pfe:$PLEIADESDIR/figs_R0/$figs ./$LOCALDIR/
        done
    done
done

