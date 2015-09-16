# HPC Coursework

Base coursework for HPC 2015 class

* Source code for simulation is in the simulation/ folder
* Inputs are in the inputs/ folder
* Results checking scripts are in the check/ folder

## Compiling and running

To compile, go into the simulation folder and type `make`. Editing the values for `CC` and `CFLAGS` in the Makefile can be used to enable different compiler options or use a different compiler. These can also be passed on the command line:

    $ make CFLAGS="-O3 -fopenmp"

Inputs, final state output files, and parameter files are all specified on the command line of the `lbm` executable. eg:

    $ ./lbm -a av_vels.dat -f final_state.dat -p ../inputs/input_128x128.params

The options for the program can be found by passing the `-h` flag to the executable.

It should be possible to compile and run the program in one by typing `make run `. By default there are some standard settings in the makefile which control which input files to use and where to output to. To use different input files, this can also be specified on the `make` command line like such:

    $ make run FINAL_STATE_FILE=./other_output_file.dat PARAM_FILE=./my_input_file.params

By default, all of these options will use the 128x128 input parameter file and the associated reference results.

## Checking results

Running `make check` will check the output file (average veclocities and final state) against some reference results. By default, it should look something like this:

    $ make check
    python ../check/check.py --ref-av-vels-file=../check/input_128x128.av_vels_ref.dat --ref-final-state-file=../check/input_128x128.final_state_ref.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
    Total difference in av_vels : 5.270812566515E-11
    Biggest difference (at step 1219) : 1.000241556248E-14
      1.595203170657E-02 vs. 1.595203170658E-02 = 6.3e-11%

    Total difference in final_state : 5.962977334129E-11
    Biggest difference (at coord (6,2)) : 1.000588500943E-14
      3.329122639178E-02 vs. 3.329122639179E-02 = 3e-11%

    Both tests passed!

This script takes both the reference results and the results to check (both average velocities and final state). This is also specified in the makefile and can be changed like the other options:

    $ make run PARAM_FILE=./my_input_file.params
    ./lbm -a ./av_vels.dat -f ./final_state.dat -p ../inputs/input_128x128.params
    ...
    $ make check REF_FINAL_STATE_FILE=./my_final_state_ref.dat REF_AV_VELS_FILE=./my_av_vels_ref.dat
    python ../check/check.py --ref-av-vels-file=./my_av_vels_ref.dat --ref-final-state-file=./my_final_state_ref.dat --av-vels-file=./av_vels.dat --final-state-file=./final_state.dat
    ...

All the options for this script can be examined by passing the --help flag to it.

    $ python ../check/check.py --help
    usage: check.py [-h] [--tolerance TOLERANCE] --ref-av-vels-file
                    REF_AV_VELS_FILE --ref-final-state-file REF_FINAL_STATE_FILE
    ...

A plot of the final state can also be done by doing `make plot`, which will create a .png file in the directory. _NOTE: This only works with a final state file called final\_state.plt_. If you want to run this on your own computer, you need gnuplot installed.

