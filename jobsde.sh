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
../abd2d --rei=0 --niu-type=0 --maxwell-grid=50 --find-grid=24 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..
# No. 6
runpath=m25n40de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=25 --find-grid=40 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..
# No. 7
runpath=m75n13de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=75 --find-grid=13 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..
# No. 8
runpath=m100n10de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=100 --find-grid=10 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..
# No. 9
runpath=m50n25de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=50 --find-grid=25 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..
# No. 1
runpath=m50n2de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=50 --find-grid=2 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..

# No. 2
runpath=m50n6de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=50 --find-grid=6 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..
# No. 3
runpath=m50n10de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=50 --find-grid=10 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..
# No. 4
runpath=m50n20de
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=50 --find-grid=20 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee co.txt 2>&1 > /dev/null &
cd ..
