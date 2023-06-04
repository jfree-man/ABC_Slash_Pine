#!/bin/bash
# Here you should provide the sbatch arguments to be used in all jobs in this serial farm
# It has to contain the runtime switch (either -t or --time):
#SBATCH --time=00-01:30:00
#SBATCH --mem=500M
#  You have to replace Your_account_name below with the name of your account:
#SBATCH --account=def-bolker
#SBATCH --mail-type=ALL
#SBATCH --mail-user=freemjc@mcmaster.ca

# Don't change this line:
task.run
