/* 
 * File:   abd2d.c
 * Author: skiloop
 *
 * Created on 2013年5月9日, 下午4:44
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "initial.h"
#include "fdtd.h"
#include "dataType.h"
#include "commonData.h"
#include "pml.h"
#include "connectingInterface.h"

void input(int argc, char*argv[]);

/*
 * 
 */
int main(int argc, char** argv) {

    input(argc, argv);
    initCommonData();
    CalDomainSize();
    InitBndCtrl();
    InitCoeff();
    InitPMLBound();
    //InitPMLBoundi();
    IntTtlFldDmnBnd();
    //InitIncFdtd();
    fdtd();
    FreeDelayArrays();
    FreePMLSpace();
    SaveToFile();
    FreeSpace();
    //FreePMLSpace();
    return (EXIT_SUCCESS);
}

#define MAX_BUFFER 201

void input(int argc, char*argv[]) {
    void PrintInput();
    void help();
    int i;
    //int cnt;
    //char buffer[MAX_BUFFER], *pstr;
    //MyDataF tdoub;
    if (argc < 2) {
        fprintf(stderr, "Invalid input!\n");
        help();
        exit(-1);
    }
    for (i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--niu-type=", 11) == 0) {
            niutype = atoi(argv[i] + 11);
        } else if (strncmp(argv[i], "--matlab=1", 10) == 0) {
            IsMatlabSim = 1;
        } else if (strncmp(argv[i], "--rei=", 6) == 0) {
            rei = atof(argv[i] + 6);
        } else if (strncmp(argv[i], "--niu-e=", 8) == 0) {
            mu_e = atof(argv[i] + 8);
        } else if (strncmp(argv[i], "--tm=1", 6) == 0) {
            IsTMx = 0;
            IsTEx = 1;
        } else if (strncmp(argv[i], "--maxwell-grid=", 15) == 0) {
            maxwellGridSize = atoi(argv[i] + 12);
            if (maxwellGridSize < 4)maxwellGridSize = MAXWELL_MESH_SIZE;
        } else if (strncmp(argv[i], "--fine-grid=", 12) == 0) {
            m = atoi(argv[i] + 12);
            if (m < 4)m = FINE_GRID_SIZE;
        } else if (strncmp(argv[i], "--total-time=", 13) == 0) {
            totaltime = atof(argv[i] + 13);
            if (totaltime <= 0)totaltime = TOTAL_TIME;
        } else if (strncmp(argv[i], "--is-connect=", 13) == 0) {
            isConnect = atoi(argv[i] + 13);
        } else if (strncmp(argv[i], "--with-density=", 15) == 0) {
            IfWithDensity = atoi(argv[i] + 15);
        } else if (strncmp(argv[i], "--e-max=", 8) == 0) {
            E_0 = atof(argv[i] + 8);
        }
#ifdef _OPENMP
        else if (strncmp(argv[i], "--thread-count=", 15) == 0) {
            thread_count = atoi(argv[i] + 15);
        }
#endif
    }
    if (niutype == 4) {
        if_erms_E_max = 1;
    }
    PrintInput();
}

void help() {
    printf("Usage:\n");
    printf("XXXX option\n\n");
    printf("--niu-type=n\tChoose niu formula type to compute.\n");
    printf("\t1\t------\tMorrow and Lowke's\n");
    printf("\t2\t------\tNikonov's\n");
    printf("\t3\t------\tKang's\n");
    printf("\tothers\t------\tdefault\n");
    printf("--matlab=[1,0]\tUse Matlab engine to show results\n");
    printf("--rei=XXX\tSet rei\n");
    printf("--niu-e=XXX\tSet Niu_e\n");
    printf("--tm=1\tset to TM mode\n");
    printf("--total-time=[>0]\ttotal simulation time\n");
    printf("--fine-grid=\thow many fine grid cells per Maxwell cell\n");
    printf("--maxwell-grid=\thow many Maxwell cells per wavelength\n");
    printf("--is-connect=[0,1]\tuse connecting interface or not\n");
    printf("--with-density=[0,1]\twether with density\n");
    printf("--e-max= \tset E field amptidute\n");
#ifdef _OPENMP
    printf("--thread-count=n\tset number of threads to run the job\n");
#endif
}

void PrintInput() {
    printf("--niu-type=%d\n", niutype);
    printf("--matlab=%d\n", IsMatlabSim);
    printf("--rei=%5.3e\n", rei);
    printf("--niu-e=%5.3e\n", mu_e);
    printf("--tm=%d\n", IsTMx);
    printf("--fine-grid=%d\n", m);
    printf("--total-time=%5.3e\n", totaltime);
    printf("--is-connect=%d\n", isConnect);
    printf("--with-density=%d\n", IfWithDensity);
    printf("--e-max=%5.3e\n", E_0);
#ifdef _OPENMP
    printf("--thead-count=%d\n", thread_count);
#endif
}

