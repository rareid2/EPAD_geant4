#!/bin/bash
#SBATCH --nodes=1                            # Number of requested nodes
#SBATCH --time=1-10:59:00                   # Max wall time
#SBATCH --qos=blanca-lair                    # Specify testing QOS
#SBATCH --ntasks=40                          # Number of tasks per job
#SBATCH --job-name=Geant4_EPP_INV            # Job submission name
#SBATCH --output=GEANT4-EPP-INV.%j.out   # Output file name with Job ID


# Written by:	Grant Berland
# Date:		28 December 2018
# Purpose: 	This script submits the Geant4 shielding simulation N times to the Slurm job scheduler


# load any modules needed to run your program  
ml gcc/10.2.0

# Re-source the build paths
source /projects/$USER/gbuild/bin/geant4.sh

# Run your program
cd /projects/$USER/EPAD_geant4/build/

./main ../macros/run_parallel_beam.mac