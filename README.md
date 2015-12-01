inabd2d
=======

inabd2d

#REQUIRE
C compiler with openmp

#HOW TO RUN ON LINUX
* 1 run make to compile
* 2 run run.sh scripts for simple running, you can change the parameters defined in run.sh to use different parameters for EM and FDTD condition

#HOW TO RUN ON WINDOWS
* 1 compile with your windows compiler
* 2 run the executable file 

#PARAMETES TO RUN  THE PROGRAM

* --rei> > > > > > recombination coefficient
* --nui-type> > > way to calculate ironlization coefficients
* --maxwell-grid> how many Maxwell grid per wavelength
* --fine-grid> > > how many density grid per Maxwell grid
* --tm> > > > > > wave type 1 for TM wave, others for TE wave
* --with-density> 1 for non-free space, others for free space
* --is-connect> > 1 for using connecting interface, others not
* --total-time> > how long the the program simulates in T


