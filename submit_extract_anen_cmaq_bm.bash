#!/bin/bash -l

#PBS -N extract_anen_cmaq_bm
#PBS -A NWSA0002
#PBS -l select=1:ncpus=1:mem=10GB
#PBS -l walltime=10:00
#PBS -q casper
#PBS -j oe

module load conda/latest
conda activate npl-2022b

python extract_anen_cmaq_bm_at_airnow_use_val_leadtime.py 20210107_06
