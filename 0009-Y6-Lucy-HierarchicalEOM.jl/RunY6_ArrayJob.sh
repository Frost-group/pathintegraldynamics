#!/bin/bash
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=32:mem=250gb
#PBS -J 1-6
# Load modules for any applications

module load Julia

# Commands for program I want to run

cd $HOME/HEOM/Y6
# First argument is dimer number
# Second argument is tier of the hierarchy to stop at 
julia -t 32 RunY6HEOM.jl $PBS_ARRAY_INDEX 3