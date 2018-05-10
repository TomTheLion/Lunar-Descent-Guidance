# Lunar-Descent-Guidance

This readme is a work in progress, for now I just wanted to have a quick reference.

The KOS files are in the LDG folder
    The ldg_func.ks file contains the main functions of the algorithm.
    The ldg_math.ks file contains general math and matrix manipulation functions.
    The ldg_cser.ks file is a conic state extrapolation routine that I borrowed from Noiredd's PEGAS program, this is pand5461's implementation of the conic state extrapolation.
    The ldg.ks file loads and runs the rest of the program.

A python simulation is in the Simulation folder. This simulation is required to determine the input variables of the Lunar Descent Algorithm specifically RBRIGZ, JBRFGZ, SBRFGX, and SBRFGZ.
