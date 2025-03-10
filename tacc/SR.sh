#!/bin/bash

k=$(($1 + 1))

echo $k
Rscript SR-simulation.R -cores 1 -reps 2000 -batch $k
