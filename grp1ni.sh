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

# No. 8
runpath=m100n10ni
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=2 --maxwell-grid=100 --find-grid=10 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..

# No. 1
runpath=m50n2ni
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=2 --maxwell-grid=50 --find-grid=2 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..

# No. 2
runpath=m50n6ni
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=2 --maxwell-grid=50 --find-grid=6 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..
