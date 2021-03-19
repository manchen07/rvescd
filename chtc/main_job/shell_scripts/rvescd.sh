#!/bin/bash

k=$(($1 + 1))

echo $k
Rscript RVE-simulation.R -cores 1 -reps 2400 -batch $k
