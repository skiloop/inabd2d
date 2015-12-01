#! /bin/sh

# jobs:
# const dxM
# No.	n	m
# 1	50	2
# 2	50	6
# 3	50	10
# 4	50	20
# 5	50	24
# const 
# No.	n	m
# 6	25	40
# 5	50	20
# 7	75	13
# 8	100	10
# compare:
# 9	50	25

# ===========================
#           runjobs 
# ===========================

# No. 5
runpath=m50n24de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=50 --fine-grid=24 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..

# No. 6
runpath=m25n40de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=25 --fine-grid=40 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..

