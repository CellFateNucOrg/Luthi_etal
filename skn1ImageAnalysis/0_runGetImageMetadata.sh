#!/bin/bash
#SBATCH --job-name="score images"
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=1

source $CONDA_ACTIVATE devbio-napari-env


work_dir="/mnt/external.data/MeisterLab/pmeister/skn-1_deletions/240411_skn-1_deletions"
#work_dir='/Users/semple/Documents/MeisterLab/otherPeopleProjects/fountains/microscopy/20240411_skn-1_deletions'

./getImageMetadata.py -d $work_dir



