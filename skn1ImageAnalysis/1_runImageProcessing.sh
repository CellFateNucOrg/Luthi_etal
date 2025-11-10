#!/bin/bash
#SBATCH --job-name="score images"
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=1

source $CONDA_ACTIVATE devbio-napari-env


python ./project_getStats_range1.py



