for d in work*/
do
    cd $d
    png2mp4.sh fluid_parameter_plots/ fluid_parameter_plots.mp4 30
    png2mp4.sh figs_R0/R0_ figs_R0.mp4 30
    cd $OLDPWD
done
