
/*
 *      This file defines functions that update scatter fields.
 *
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#ifdef _OPENMP
#include<omp.h>
#endif
#include "common.h"
#include "initial.h"
#include "fdtd.h"
#include "density.h"
#include "dataType.h"
#include "dataSaving.h"
#include "pml.h"
#include "commonData.h"
#include "connectingInterface.h"
#include "matlabSimulation.h"


///////////////////////////////
// FDTD DATA
///////////////////////////////
//EM field
MyStruct Ey, Ex, Ez, Hx, Hy, Hz;
MyStruct Ey_pre, Ex_pre, Ez_pre;
MyStruct Vey, Vex, Vez;
MyStruct ne, ne_pre;
MyStruct beta;
MyStruct Erms;
MyStruct Ermsx, Ermsy;
MyStruct Nu_c; // 循环利用的碰撞率

//Coeffients
MyDataF Chxez, Chyez, Chzex, Chzey;
MyStruct Cevx, Cevy, Cevz;
MyStruct Ceex, Ceey, Ceez;
MyStruct Cehx, Cehy, Cehz;

MyDataF Cve;
MyStruct Cvvx, Cvvy, Cvvz;
MyStruct Cvex, Cvey, Cvez;

///////////////////////////////
//formation 3
///////////////////////////////

void UpdateEField() {

    int i, j, index;
    int indx, indy;
    //MyDataF Eztmp;
    if (IsTMx) {
        BackupMyStruct(Ex, Ex_pre);
        BackupMyStruct(Ey, Ey_pre);
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,indx)
#endif
        for (i = 0; i < Ex.nx; i++) {
            for (j = pjs; j < pje; j++) {
                index = i * Hz.ny + j;
                indx = i * Ex.ny + j;
                Ex.data[indx] = Ceex.data[indx]*(Ex.data[indx]) +
                        Cevx.data[indx] * Vex.data[indx] + Cehx.data[indx]*
                        (Hz.data[index] - Hz.data[index - 1]);
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,indy)
#endif
        for (i = pis; i < pie; i++) {
            for (j = 0; j < Ey.ny; j++) {
                indy = i * Ey.ny + j;
                Ey.data[indy] = Ceey.data[indy]*(Ey.data[indy])
                        + Cevy.data[indy] * Vey.data[indy] + Cehy.data[indy]*
                        (Hz.data[indy] - Hz.data[indy - Hz.ny]);
            }
        }
        //AdjustEFieldAtCnntIntfc();
    }
    if (IsTEx) {
        BackupMyStruct(Ez, Ez_pre);
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,indx,indy)
#endif
        for (i = pis + 1; i < pie; i++) {
            for (j = pjs + 1; j < pje; j++) {
                index = i * Ez.ny + j;
                indy = i * Hy.ny + j;
                indx = i * Hx.ny + j;

                Ez.data[index] = Ceez.data[index]*(Ez.data[index]) +
                        Cevz.data[index] * Vez.data[index] + Cehz.data[index]*(
                        (Hy.data[indy] - Hy.data[indy - Hy.ny])-
                        (Hx.data[indx] - Hx.data[indx - 1]));

            }
        }
    }
}

void UpdateMField() {

    int i, j, index, ind;
    //MyDataF tmp1,tmp2,tmp3;
    if (IsTEx) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,ind)
#endif
        for (i = 0; i < Hx.nx; i++) {
            for (j = pjs; j < pje; j++) {
                index = i * Hx.ny + j;
                ind = i * Ez.ny + j;
                Hx.data[index] = Hx.data[index] + Chxez * (Ez.data[ind + 1] - Ez.data[ind]);
            }
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,ind)
#endif
        for (i = pis; i < pie; i++) {
            for (j = 0; j < Hy.ny; j++) {
                index = i * Hy.ny + j;
                ind = i * Ez.ny + j;
                Hy.data[index] = Hy.data[index] + Chyez * (Ez.data[ind + Ez.ny] - Ez.data[ind]);
            }
        }
    }
    if (IsTMx) {
        int ind2;
#ifdef _OPENMP
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,ind,ind2)
#endif
        for (i = pis; i < pie; i++) {
            for (j = pjs; j < pje; j++) {
                index = i * Hz.ny + j;
                ind = i * Ex.ny + j;
                ind2 = i * Ey.ny + j;
                Hz.data[index] +=
                        Chzex * (Ex.data[ind + 1] - Ex.data[ind]) +
                        Chzey * (Ey.data[ind2 + Ey.ny] - Ey.data[ind2]);
#ifdef _DEBUG
                if (index == Hz.nx * Hz.ny / 2 + Hz.ny / 2)
                    j = j;
#endif
            }
        }
    }
}


int Step = 1;
int MTimeStep = 1;
int cnepx, cnepy;

void AdBound();
void InitSources();

void fdtd() {

    int CurTimeStep; //current time step
    int xpos, ypos, sxpos, sypos;
    //FILE *fne;
    //int i,j;
    //MyDataF *CapEF;
    void PrintSourceSize(int, int, int, int);
    MyDataF cdtF;
    //int MTimeStep=1,MultiSize;
    //int cnt=0;
    int step_per_half_ns8 = (int) (0.5 + 0.125e-9 / dt_F);
    clock_t t_start, t_end;
    //MyDataF RealEz,rhtb;
    int stpos = NUMBER_OF_CELLS_IN_PML_BOUND + SCATTER_FIELD_DOMAIN_BND_SIZE;
    int index;
    sxpos = (int) (0.5 + stpos + (int) (2.25 / 3.0 * (nx - 2 * stpos))); //x position of sources
    sypos = (int) (0.5 + ny / 2);
    if (IsTEx) {
        index = sxpos * Ez.ny + sypos;
    } else {
        index = sxpos * Hz.ny + sypos;
    }
    //fne=fopen("cnep.dat","w");
    //y position of sources
    cnepx = m * ((int) (0.5 + 2.25 * lamda / dx) + tpis);
    cnepy = m * (nyp1 / 2);
    xpos = (int) (0.5 + nx / 2 + lamda / dx); //x position of field to be captured
    ypos = (int) (0.5 + ny / 2 + lamda / dy); //y position of field to be captured

    PrintSourceSize(xpos, ypos, sxpos, sypos);
    MultiSize = (int) (0.5 + dt_F / dt);

    printf("\ncnepx = %d\tcnepy = %d\n\n", cnepx, cnepy);
    printf("\nm\t=\t%d\nMultiSize\t=\t%d", m, MultiSize);
    printf("\nNumber of time steps: %d\n", TotalTimeStep);
    printf("Density Time Step: %d\n", Density_Time_Step * m);
    printf("Step per half ns: %d\n", step_per_half_ns8);
    printf("\n*************************************************************\n");

    sxpos = (int) (0.5 + nx / 2 + 0.125 * lamda / dx);
    sypos = (int) (0.5 + (tpjs + tpje) / 2);
    CurTime = -half_dt;

    initSource();
    InitSim();
    //InitIncFdtd();
    if (isConnect)initconnect(); ////InitFields();////Init fields before marching  loop i.e. at t = 0
    t_start = clock();
    for (CurTimeStep = 1; CurTimeStep <= Density_Time_Step && TotalTimeStep >= Step; CurTimeStep++) {

        for (cdtF = 0, MTimeStep = 1; MTimeStep <= MultiSize && cdtF < dt_F && TotalTimeStep >= Step; cdtF = cdtF + dt, MTimeStep++, Step++) {
            //MyDataF ezl,ezlr;

            CurTime += half_dt;

            UpdateMField();
            //connecting interface
            if (isConnect)mconnect(CurTime);
            //DispEMFields(CurTimeStep);
            //Hz_s.data[Hz_s.ny*Hz_s.ny/2+Hz_s.ny/2] += Source(CurTime);
            //UpdMFieldForPML();
            UpdateMFieldForPML(Hx, Hy, Hz, Ex, Ey, Ez);
            //DispEMFields(CurTimeStep);

            CurTime += half_dt;
            //DispEMFields(CurTimeStep);
            UpdateEField();
            if (!isConnect) {
                if (IsTEx)Ez.data[index] += E0 * Source(CurTime);
                if (IsTMx)Hz.data[index] += H0 * Source(CurTime);
            }
            //DispEMFields(CurTimeStep);
            //connecting interface
            if (isConnect)econnect(CurTime);
            //Adjust_E_Field(CurTime-half_dt);////AddEInc(CurTime);//
            //DispEMFields(CurTimeStep);
            //UpdEFieldForPML();
            UpdateEFieldForPML(Hx, Hy, Hz, Ex, Ey, Ez);
            //DispEMFields(CurTimeStep);
            if (IfWithDensity) {
                UpdateVelocity();
                //CaptureE(CapEF,Step,xpos,ypos);
                if (if_erms_E_max == 1) {
                    CalSumESqrt_Emax();
                } else {
                    CalSumESqrt();
                }
            }
            if (Step % 50 == 0) {
                DispEMFields(CurTimeStep);
            }
        }
        SaveCapField(Step); // 保存数据
        if (IfWithDensity) {
            if (if_erms_E_max != 1) {
                CalErmsAtCoarseGrid();
            }
            InterpolatErms();
            UpdateDensity();
            UpdateCoeff();
            //SaveErms(CurTimeStep);
            //printf("ne = %15.14e\n",ne.data[cnepx*ne.ny+cnepy]);
            DispNe(CurTimeStep);
            //DispEMFields(CurTimeStep);
            ResetStructData(Erms);
        }
        if (IsTEx) {
            printf("%d\t%5.4e\t%f ns\n", Step, Ez.data[index], CurTime / 1e-9);
        } else {
            printf("%d\t%5.4e\t%f ns\n", Step, Hz.data[index], CurTime / 1e-9);
        }
        ////
        //if(CurTimeStep%1==0)
    }//END FOR
    //SaveD(cnt);
    t_end = clock();
    printf("Total time used : %ld\n", t_end - t_start);
    //system("pause");
    EndMatlabSim();
    ////Save Sample Field
    //scpfld(CapEF,TotalTimeStep);
    //free(CapEF);
    end_connect();
    free(srcdat);

}//END OF FDTDLOOP

void PrintSourceSize(int xpos, int ypos, int sxpos, int sypos) {
    printf("\nnx\t=\t%d\nny\t=\t%d\nxpos\t=\t%d\nypos\t=\t%d\nsxpos\t=\t%d\nsypos\t=\t%d\nnbound\t=\t%d\n", nx, ny, xpos, ypos, sxpos, sypos, nbound);
}

void FreeSpace() {
    if (IsTMx) {
        freeData(&Ex);
        freeData(&Ex_pre);

        freeData(&Vex);
        freeData(&Cevx);
        freeData(&Ceex);
        freeData(&Cehx);

        freeData(&Ey);
        freeData(&Ey_pre);

        freeData(&Vey);
        freeData(&Cevy);
        freeData(&Ceey);
        freeData(&Cehy);
        freeData(&Hz);
    }
    if (IsTEx) {
        freeData(&Hx);
        freeData(&Hy);
        freeData(&Ez);
        freeData(&Ez_pre);

        freeData(&Vez);
        freeData(&Ceez);
        freeData(&Cevz);
        freeData(&Cehz);
    }
    freeData(&Cvvx);
    freeData(&Cvvy);
    freeData(&Cvvz);
    freeData(&Cvex);
    freeData(&Cvey);
    freeData(&Cvez);
    freeData(&ne);
    freeData(&Erms);
    freeData(&ne_pre);
    freeData(&beta);
    freeData(&Nu_c);
}

