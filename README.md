# Files for "Observed Extra Mixing Trends in Red Giants are Reproduced by the Reduced Density Ratio in Thermohaline Zones"

This repository contains most of the files used to create the content in the paper "Observed Extra Mixing Trends in Red Giants are Reproduced by the Reduced Density Ratio in Thermohaline Zones" by A. Fraser, M. Joyce, E. Anders, J. Tayar and M. Cantiello.


## Running the MESA grids

The code (inlists, etc.) used to run the various MESA grids in this paper can be found in the mesa_models/ folder.
The workflow for creating the main grids of models (spanning 4 thermohaline mixing prescriptions in a [Fe/H] x Mass grid) is as follows:

1. `cd mesa_models`, then run `python3 generate_run_directories.py`. This creates a new folder ("work_thermohaline_...") for each choice of mixing prescription, mass, metallicity, etc. To do so, a copy of the `template_work_thermohaline` folder is created, and the MESA inlist is modified with the appropriate choices as specified in the Python file.
2. For each run directory `WORK_DIR`, cd `WORK_DIR` then follow the full workflow laid out in the `evan_pleiades_job_submit` file. This includes:
    1.  Cleaning, making, and running the MESA model (`./clean`, `./mk`, `./rn`). Sometimes things go wrong or weird with the copying or syncing to github, so you may have to make these files executable (e.g., `chmod +x mk`). We ran all models using MESA version `mesa-r21.12.1`.
    2.  Decreasing the amount of disk space required to store the MESA outputs (`python3 mesa_to_hdf5.py`, `python3 hdf5_to_lite_hdf5.py`). The first command creates an hdf5 file of all of the MESA outputs, which reduces disk space of outputs by a factor of ~5 (filename `hdf5_LOGS.h5`). The second file further reduces disk space by only preserving some of the profile fields for only a small portion of the star around the relevant region for the analysis in this paper (`hdf5_lite.h5`).
    3.  Post-processing the lite hdf5 outputs (`python3 find_r_vs_time.py`), which produces the movie shown in e.g., figure 7 of the paper as well as a file `r_vs_time.h5` which contains time traces of r, R0, tau, log g, etc.
    4.  Removal of the original MESA LOGS/ directory to clean up disk space.
    5.  Creation of movies as in e.g., fig 7 using a the png2mp4 script originally created by Keaton Burns.
3. After this process is completed and `r_vs_time.h5` is created, we move the work directory into the post directory (`cd ..`, `cd WORK_DIR post/`).
4. Once the full grid has been moved into post/, we run `python3 generate_csv_data.py`, which generates the various csv files which are stored in the `Empirical-Magnetic-Thermohaline/figures/mesa_spread.csv` folder.

To conduct the resolution test shown in appendix B, all of the same steps were followed, but using the `generate_run_directories.py` file in `mesa_models/resolution_test/`. In step 3, finished work directories are moved to `post/resolution_test/`, and the csv file is created using the `generate_csv_data.py` file there. The output csv files are moved to `figures/resolution_test/`.

The python workflow depends on the `numpy`, `scipy`, and `hdf5` packages.

## Creating figures

To create figures 1, 2, 4-6 in the manuscript, navigate to the `figures/` directory then run `./generate_figures.sh`. 
Figure 3 is created using elder magic (IDL) and files to create it were uploaded by J. Tayar in commit [3736ee6](https://github.com/afraser3/Empirical-Magnetic-Thermohaline/commit/3736ee64cef2af9f285e81acc35ff172217fa902).
Figure 7 is created during the MESA post-processing pipeline described above.
