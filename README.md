# Value of Information for Social Networks in Conservation Planning

This repository contains the code and data for the manuscript "Fundamental insights on when social network data is most critical for conservation planning". The code to run the optimisation is dependent on Linux, while the evaluation and statistical analysis will run on Windows or Linux (and probably also iOS, but this has not been tested).   

See: [![DOI](https://zenodo.org/badge/147589655.svg)](https://zenodo.org/badge/latestdoi/147589655)

## Instructions

1) Download the code into a suitable location, maintaining the file structure in the repository

2) Install required dependencies (see below)

3) Run the optimisation for the simulations and for the fishery by setting the working directory in Matlab to "/matlab_code" and executing "run.m". Not that you need to run Matlab from within the Symbolic Perseus root directory (see notes on dependencies below).

4) Run the statistical analysis for the simulations and for the fishery by setting the working directory in R to "/r_code" and executing "analysis_CPSN_CP.r" for the simulations or "analysis_CPSN_CP_fishery.r" for the fishery. Note that since the summary outputs used in the manuscript are currently provided in the "/problems/ConsMDP/results" directory this analysis can be run without having to re-run the optimisation and evaluation (step 3 above). However, running the optimisation in step 3 above will over-write the existing summary results in this directory.     

5) Results of the optimisation and evaluation will be saved to the "/problems/ConsMDP" directory.

6) Summarised outputs from the optimisation and evaluation used in the manuscript can be found in the "/problems/ConsMDP/results" directory.     

## Dependencies

The code has the following dependencies:

1) Matlab

2) Symbolic Perseus (https://cs.uwaterloo.ca/~ppoupart/software.html). Note that when running Matlab, run it from within the Symbolic Perseus root directory once compiled. The source code for Symbolic Perseus is also provided in this repository in the "/dependencies" directory.

3) SPUDD 3.6.2 (https://cs.uwaterloo.ca/~jhoey/research/spudd/index.php). Note that prior to compiling, change 'sprintf(temp,"%s",lnames[(int) (lval-1)]);' on line 98 in the "dumpdot.c" file in the "aconvert" function to 'sprintf(temp,"%s",lnames[(int) (lval-1+0.5)]);'. This fixes a rounding error. The source code for SPUDD is also provided in this repository in the "/dependencies" directory. You may also need to edit some other aspects of the code for SPUDD to be able to compile it depending on the version of Linux and the complier used.      

4) R (https://www.r-project.org/)

5) Jags (http://mcmc-jags.sourceforge.net/)
