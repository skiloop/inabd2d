#! /bin/sh

# ===========================
#           runjobs 
# ===========================

# No. 1
runpath=nikonov
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=2 --maxwell-grid=50 --find-grid=10 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee r.log 2>&1 > /dev/null 
cd ..
#
# No. 2
runpath=morrow
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=1 --maxwell-grid=50 --find-grid=10 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee r.log 2>&1 > /dev/null 
cd ..

 No. 3
runpath=kang
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=3 --maxwell-grid=50 --find-grid=10 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee r.log 2>&1 > /dev/null 
cd ..

# No. 4
runpath=default
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=0 --maxwell-grid=50 --find-grid=10 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee r.log 2>&1 > /dev/null 
cd ..

# No. 5
runpath=formula4
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=4 --maxwell-grid=50 --find-grid=10 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee r.log 2>&1 > /dev/null 
cd ..

# No. 6
runpath=formula5
mkdir -p $runpath
cd $runpath
../abd2d --rei=0 --niu-type=5 --maxwell-grid=50 --find-grid=10 --is-connect=1 --with-density=1 --tm=1 --total-time=150 | tee r.log 2>&1 > /dev/null 
cd ..

