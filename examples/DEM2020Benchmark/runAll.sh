#!/bin/bash

# This script will download and run all cases from the 2020 benchmark

# PREREQUISITES:
# 1. an internet connexion (else please comment out the download part and provide the scripts in current path)
# 2. a recent version of yade (not earlier than 01/2021 else one python helper for converting mesh files will be missing)

# The scripts themselves 
# Cases with output provided in 02/2021 are signaled by ending the line [*]


# That url should point to valid gitlab branch
export YADE_BRANCH="https://gitlab.com/yade-dev/trunk/-/raw/edea1b5e219691753662104bb7cfac0986536aed/examples/DEM2020Benchmark"
# latest would be:
# export YADE_BRANCH=https://gitlab.com/yade-dev/trunk/-/raw/master/examples//DEM2020Benchmark

wget $YADE_BRANCH/Case1_SiloFlow.py
wget $YADE_BRANCH/Case2_rotating_drum.py
wget $YADE_BRANCH/Case3_PenetrationTest.py

export YADE="yadedaily"
# export YADE="/path/to/my/own/yadeVersion"

export OMP_THREADS=4 # OpenMP threads, should be less than number of cores. 12 max suggested.

$YADE -j $OMP_THREADS -n -x Case1_SiloFlow.py small M1 #[*]
$YADE -j $OMP_THREADS -n -x Case1_SiloFlow.py small M2
$YADE -j $OMP_THREADS -n -x Case1_SiloFlow.py large M1 #[*]
$YADE -j $OMP_THREADS -n -x Case1_SiloFlow.py large M2
$YADE -j $OMP_THREADS -n -x /mnt/Case2_rotating_drum.py #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest_SI.py 25000 M1 #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest_SI.py 50000 M1 #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest_SI.py 100000 M1 #[*]
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest_SI.py 25000 M2
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest_SI.py 50000 M2
$YADE -j $OMP_THREADS -n -x Case3_PenetrationTest_SI.py 100000 M2

$YADE plotBenchmark.py # would run just as well with python3
