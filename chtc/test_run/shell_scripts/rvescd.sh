#!/bin/bash

# k starts from 0 in condor, make it start from 1
k=$(($1 + 1))

echo $k
Rscript RVE-simulation.R -cores 1 -reps 2400 -batch $k
