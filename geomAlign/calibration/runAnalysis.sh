#!/bin/bash  
echo "Running gmake"
cd /star/u/jdb/WORK/TOFCalib/doCalib
cd ana
gmake

starver SL14a

echo "Running analysis on ntuples"
./analysis /star/u/jdb/WORK/TOFCalib/run13Align/output
#mv align.root ../fit/