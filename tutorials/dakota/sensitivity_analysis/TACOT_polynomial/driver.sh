#!/bin/bash

# $1 and $2 are special variables in bash that contain the 1st and 2nd 
# command line arguments to the script, which are the names of the
# Dakota parameters and results files, respectively.

params=$1
results=$2

############################################################################### 
##
## Pre-processing Phase -- Generate/configure an input file for your simulation 
##  by substiting in parameter values from the Dakota paramters file.
##
###############################################################################

dprepro $params constantProperties.template constantProperties.i

############################################################################### 
##
## Execution Phase -- Run your simulation
##
###############################################################################


./Allclean_pato > /dev/null
cp constantProperties.i TACOT-UQ/constantProperties
 ./Allrun_pato > /dev/null 
rm constantProperties.i

############################################################################### 
##
## Post-processing Phase -- Extract (or calculate) quantities of interest
##  from your simulation's output and write them to a properly-formatted
##  Dakota results file.
##
###############################################################################

cat results > $results 

