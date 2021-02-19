#!/bin/bash

# This script will download and run all cases from the 2020 benchmark

# PREREQUISITES:
# 1. an internet connexion (else please comment out the download part and provide the scripts in current path)
# 2. a recent version of yade (not earlier than 01/2021 else one python helper for converting mesh files will be missing)

# The scripts themselves 
# Cases with output provided in 02/2021 are signaled by ending the line [*]


# That url should point to valid gitlab branch
export YADE_BRANCH="https://gitlab.com/yade-dev/trunk/-/raw/99f912067159bf8a475ea217c8c46db94f67215d/examples/DEM2020Benchmark"
# latest would be:
# export YADE_BRANCH=https://gitlab.com/yade-dev/trunk/-/raw/master/examples//DEM2020Benchmark


wget $YADE_BRANCH/Case1_SiloFlow.py
wget $YADE_BRANCH/Case2_rotating_drum.py
wget $YADE_BRANCH/Case3_PenetrationTest.py

yadedaily -j 4 -n -x Case1_SiloFlow.py small M1 #[*]
yadedaily -j 4 -n -x Case1_SiloFlow.py small M2
yadedaily -j 4 -n -x Case1_SiloFlow.py large M1 #[*]
yadedaily -j 4 -n -x Case1_SiloFlow.py large M2
yadedaily -j 4 -n -x /mnt/Case2_rotating_drum.py #[*]
yadedaily -j 4 -n -x Case3_PenetrationTest_SI.py 25000 M1 #[*]
yadedaily -j 4 -n -x Case3_PenetrationTest_SI.py 50000 M1 #[*]
yadedaily -j 4 -n -x Case3_PenetrationTest_SI.py 100000 M1 #[*]
yadedaily -j 4 -n -x Case3_PenetrationTest_SI.py 25000 M2
yadedaily -j 4 -n -x Case3_PenetrationTest_SI.py 50000 M2
yadedaily -j 4 -n -x Case3_PenetrationTest_SI.py 100000 M2

yadedaily plotBenchmark.py
