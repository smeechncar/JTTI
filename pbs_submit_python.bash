#!/bin/bash

#PBS -N stats_by_epa_reg
#PBS -A NWSA0002
#PBS -l walltime=00:15:00
#PBS -q casper
#PBS -l select=1:ncpus=1
#PBS -j oe
#PBS -M jaredlee@ucar.edu

python stats_by_epa_region.py
