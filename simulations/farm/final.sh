#!/bin/bash

#SBATCH --time=00-01:00:00
#SBATCH --mem=4G
#SBATCH --mail-user=freemjc@mcmaster.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-bolker

module load r/4.1.2
Rscript ~/ABC_Slash_Pine/simulations/combine.R
