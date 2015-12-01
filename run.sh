#! /bin/sh

# ===========================
#           runjobs 
# ===========================


runpath=out
run=abd2d
rei=0
nuiType=2
maxwellGrid=50
connect=1
density=1
tmWave=1
runtime=150
if [ ! -f ${run} ];then
	echo "${run} does not exist,please execute 'make' to compile";
	exit;
fi
mkdir -p ${runpath}
cd ${runpath}
../${run} --rei=${rei} --niu-type=${nuiType} --maxwell-grid=${maxwellGrid} --find-grid=${densityGrid} --is-connect=${connect} --with-density=${density} --tm=${tmWave} --total-time=${runtime} > r.log 2>&1 > /dev/null &
echo "jobs is running in background"
echo "data output is in directory ${runpath}";
echo "to see if still running please execute: ps | grep ${run}";
cd ..
