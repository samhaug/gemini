#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH -o GEMINI.out
#SBATCH -e GEMINI.err
#SBATCH -p DTM
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
ulimit -s unlimited
./do_all fuck  >& OUTPUT
