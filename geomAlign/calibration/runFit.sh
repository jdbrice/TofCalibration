#!/bin/bash  

echo "Runing the Fit Scripts"
cd /star/u/jdb/WORK/TOFCalib/doCalib
cd fit
root -b -q getX0.C
root -b -q getY0.C
root -b -q genOffset.C