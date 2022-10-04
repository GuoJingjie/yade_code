#!/bin/bash

# This script will download and run all cases from the 2020 benchmark. It is standalone and should run on every linux system with yade installed.
# The script itself can retrieved here:
#   $ wget https://gitlab.com/yade-dev/trunk/-/raw/master/examples/DEM2020Benchmark/runAll.sh

# PREREQUISITES:
# 1. an internet connexion (else please comment out the download part and provide the scripts in current path)
# 2. a recent version of yade (not earlier than 01/2021 else one python helper for importing mesh files would be missing)
# installation from deb packages ('focal' is ubuntu20.04, replace with relevant name of your ubuntu/debian distro [1]):
#  
#     sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ focal main" >> /etc/apt/sources.list'
#     wget -O - http://www.yade-dem.org/packages/yadedev_pub.gpg | sudo apt-key add -
#     sudo apt-get update
#     sudo apt-get install yadedaily
#
#  [1] if needed, more install instructions here: https://yade-dem.org/doc/installation.html

# ALTERNATIVE INPUTS
# Each script is supposed to retrieve input data (particle and wall positions) from the internet. 
# If you prefer to use local files, please make sure the input files are like below relative to current path (where the scripts are executed), with the same names.
# Make sure the first line with column headers (if any) starts with a '#', and that the particles are defined with columns 'x y z rad' (in this order).

# ./inputData/
# 	Case1_SiloFlow_PartCoordinates_Case1_smallM1.txt
# 	Case1_SiloFlow_PartCoordinates_Case1_largeM1.txt
# 	Case1_SiloFlow_PartCoordinates_Case1_smallM2.txt
# 	Case1_SiloFlow_PartCoordinates_Case1_largeM2.txt
# 	Case1_SiloFlow_Walls_Case1_smallM1.txt
# 	Case1_SiloFlow_Walls_Case1_largeM1.txt
# 	Case2_Drum_Walls.txt
# 	Case2_Drum_PartCoordinates.txt
# 	Case3Walls25.txt
# 	25KParticles.txt
# 	Case3Walls50.txt
# 	50KParticles.txt
# 	Case3Walls100.txt
# 	100KParticles.txt

# OUTPUTS:
# *.png and *.txt files per-job will be in ./outputData
# synthetic png per case will be in current path, as well as a summary of the timings in timings.txt

# WARNING
# If for some reason this script execution is stopped with 'ctrl+z' before the end, then some jobs will stay behind and still use ressources, 
# don't forget 'pkill -9 yadedaily' after that (or whatever yade version it is) to really kill them


# That url should point to valid gitlab branch from where the benchmark scripts can be retrieved
export YADE_BRANCH='https://gitlab.com/yade-dev/trunk/-/raw/5cbbc95ed986b/examples/DEM2020Benchmark'
# latest would be:
# export YADE_BRANCH='https://gitlab.com/yade-dev/trunk/-/raw/master/examples//DEM2020Benchmark'

# what follows will not overwrite existing scripts if they are already in current folder
# make sure you erase them before running this script if you want fresh versions
wget -nc $YADE_BRANCH/Case1_SiloFlow.py
wget -nc $YADE_BRANCH/Case2_rotating_drum_mpi.py
wget -nc $YADE_BRANCH/Case3_PenetrationTest.py
wget -nc $YADE_BRANCH/plotBenchmark.py

export YADE='yadedaily'
# export YADE="/path/to/my/own/yadeVersion"

export OMP_THREADS=12 # OpenMP threads, should be less than number of cores. 12 max suggested.
export OMP_PROC_BIND=true # pin OMP threads to physical cores

# Official simulation times
export simulationTime1=5
export simulationTime2=5
export simulationTime3=0.01

# For testing other times
# export simulationTime1=0.00001
# export simulationTime2=0.00001
# export simulationTime3=0.00001

# Cases with output provided in 02/2021 are signaled by ending the line [*]
$YADE -j $OMP_THREADS -n -x Case1_SiloFlow.py small M1 $simulationTime1 #[*]
$YADE -j $OMP_THREADS -n -x Case1_SiloFlow.py small M2 $simulationTime1
$YADE -j $OMP_THREADS -n -x Case1_SiloFlow.py large M1 $simulationTime1 #[*]
$YADE -j $OMP_THREADS -n -x Case1_SiloFlow.py large M2 $simulationTime1
$YADE -j $OMP_THREADS -n -x Case2_rotating_drum_mpi.py $simulationTime2 #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest.py 25000 $simulationTime3 #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest.py 50000 $simulationTime3 #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest.py 100000 $simulationTime3 #[*]

$YADE -n -x plotBenchmark.py # would run just as well with python3
