#!/bin/bash

# This script will download and run all cases from the 2020 benchmark. It is standalone and should run on every linux system with yade installed.
# The script itself can retrieved here:
#   $ wget https://gitlab.com/yade-dev/trunk/-/raw/ee79223a8/examples/DEM2020Benchmark/runAll.sh

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


# OUTPUTS:
# *.png and *.txt files per-job will be in ./outputData
# synthetic png per case will be in current path, as well as a summary of the timings in timings.txt

# WARNING
# If for some reason this script execution is stopped with 'ctrl+z' before the end, then some jobs will stay behind and still use ressources, 
# don't forget 'pkill -9 yadedaily' after that (or whatever yade version it is) to really kill them


# That url should point to valid gitlab branch from where the benchmark scripts can be retrieved
export YADE_BRANCH='https://gitlab.com/yade-dev/trunk/-/raw/0402e0071/examples/DEM2020Benchmark'
# latest would be:
# export YADE_BRANCH='https://gitlab.com/yade-dev/trunk/-/raw/master/examples//DEM2020Benchmark'

# what follows will not overwrite existing scripts if they are already in current folder
# make sure you erase them before running this script if you want fresh versions
wget -nc $YADE_BRANCH/Case1_SiloFlow.py
wget -nc $YADE_BRANCH/Case2_rotating_drum.py
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
$YADE -j $OMP_THREADS -n -x Case2_rotating_drum.py $simulationTime2 #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest.py 25000 $simulationTime3 #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest.py 50000 $simulationTime3 #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest.py 100000 $simulationTime3 #[*]

$YADE -n -x plotBenchmark.py # would run just as well with python3
