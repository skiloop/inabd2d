
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
MyStruct Ey_i, Ex_i, Ez_i, Hx_i, Hy_i, Hz_i;
MyStruct Ey_s, Ex_s, Ez_s, Hx_s, Hy_s, Hz_s;
MyStruct Ey_s_pre, Ex_s_pre, Ez_s_pre;
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
        BackupMyStruct(Ex_s, Ex_s_pre);
        BackupMyStruct(Ey_s, Ey_s_pre);
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,indx)
        for (i = 0; i < Ex_s.nx; i++) {
            for (j = pjs; j < pje; j++) {
                index = i * Hz_s.ny + j;
                indx = i * Ex_s.ny + j;
                Ex_s.data[indx] = Ceex.data[indx]*(Ex_s.data[indx]) +
                        Cevx.data[indx] * Vex.data[indx] + Cehx.data[indx]*
                        (Hz_s.data[index] - Hz_s.data[index - 1]);
            }
        }
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,indy)
        for (i = pis; i < pie; i++) {
            for (j = 0; j < Ey_s.ny; j++) {
                indy = i * Ey_s.ny + j;
                Ey_s.data[indy] = Ceey.data[indy]*(Ey_s.data[indy])
                        + Cevy.data[indy] * Vey.data[indy] + Cehy.data[indy]*
                        (Hz_s.data[indy] - Hz_s.data[indy - Hz_s.ny]);
            }
        }
        //AdjustEFieldAtCnntIntfc();
    }
    if (IsTEx) {
        BackupMyStruct(Ez_s, Ez_s_pre);
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,indx,indy)
        for (i = pis + 1; i < pie; i++) {
            for (j = pjs + 1; j < pje; j++) {
                index = i * Ez_s.ny + j;
                indy = i * Hy_s.ny + j;
                indx = i * Hx_s.ny + j;

                Ez_s.data[index] = Ceez.data[index]*(Ez_s.data[index]) +
                        Cevz.data[index] * Vez.data[index] + Cehz.data[index]*(
                        (Hy_s.data[indy] - Hy_s.data[indy - Hy_s.ny])-
                        (Hx_s.data[indx] - Hx_s.data[indx - 1]));

            }
        }
    }
}

void UpdateMField() {

    int i, j, index, ind;
    //MyDataF tmp1,tmp2,tmp3;
    if (IsTEx) {
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,ind)
        for (i = 0; i < Hx_s.nx; i++) {
            for (j = pjs; j < pje; j++) {
                index = i * Hx_s.ny + j;
                ind = i * Ez_s.ny + j;
                Hx_s.data[index] = Hx_s.data[index] + Chxez * (Ez_s.data[ind + 1] - Ez_s.data[ind]);
            }
        }
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,ind)
        for (i = pis; i < pie; i++) {
            for (j = 0; j < Hy_s.ny; j++) {
                index = i * Hy_s.ny + j;
                ind = i * Ez_s.ny + j;
                Hy_s.data[index] = Hy_s.data[index] + Chyez * (Ez_s.data[ind + Ez_s.ny] - Ez_s.data[ind]);
            }
        }
    }
    if (IsTMx) {
        int ind2;
#pragma omp parallel for num_threads(thread_count) schedule(dynamic)\
	private(i,j,index,ind,ind2)
        for (i = pis; i < pie; i++) {
            for (j = pjs; j < pje; j++) {
                index = i * Hz_s.ny + j;
                ind = i * Ex_s.ny + j;
                ind2 = i * Ey_s.ny + j;
                Hz_s.data[index] +=
                        Chzex * (Ex_s.data[ind + 1] - Ex_s.data[ind]) +
                        Chzey * (Ey_s.data[ind2 + Ey_s.ny] - Ey_s.data[ind2]);
#ifdef _DEBUG
                if (index == Hz_s.nx * Hz_s.ny / 2 + Hz_s.ny / 2)
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
    sxpos = (int) (0.5 + ny / 2); //x position of sources
    sypos = (int) (0.5 + ny / 2);
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
    filedat = fopen("neiter.dat", "w");
    if (filedat == NULL) {
        fprintf(stderr, "Cannot open neiter.dat for output!\n");
        exit(EXIT_FAILURE);
    }
    CurTime = -half_dt;
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
            UpdateMFieldForPML(Hx_s, Hy_s, Hz_s, Ex_s, Ey_s, Ez_s);
            //DispEMFields(CurTimeStep);

            CurTime += half_dt;
            //DispEMFields(CurTimeStep);
            UpdateEField();
            if (!isConnect) {
                if (IsTEx)Ez_s.data[100 * Ez_s.ny + 50] += Source(CurTime);
                if (IsTMx)Hz_s.data[(pis + 10) * Hz_s.ny + 100] = Source(CurTime);
            }
            //DispEMFields(CurTimeStep);
            //connecting interface
            if (isConnect)econnect(CurTime);
            //Adjust_E_Field(CurTime-half_dt);////AddEInc(CurTime);//
            //DispEMFields(CurTimeStep);
            //UpdEFieldForPML();
            UpdateEFieldForPML(Hx_s, Hy_s, Hz_s, Ex_s, Ey_s, Ez_s);
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
            if (Step == 100) {
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
        printf("%d\t%f ns\n", Step, CurTime / 1e-9);
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
    fclose(filedat);

}//END OF FDTDLOOP

void PrintSourceSize(int xpos, int ypos, int sxpos, int sypos) {
    printf("\nnx\t=\t%d\nny\t=\t%d\nxpos\t=\t%d\nypos\t=\t%d\nsxpos\t=\t%d\nsypos\t=\t%d\nnbound\t=\t%d\n", nx, ny, xpos, ypos, sxpos, sypos, nbound);
}

void FreeSpace() {
    if (IsTMx) {
        freeData(&Ex_s);
        freeData(&Ex_s_pre);

        freeData(&Vex);
        freeData(&Cevx);
        freeData(&Ceex);
        freeData(&Cehx);

        freeData(&Ey_s);
        freeData(&Ey_s_pre);

        freeData(&Vey);
        freeData(&Cevy);
        freeData(&Ceey);
        freeData(&Cehy);
        freeData(&Hz_s);
    }
    if (IsTEx) {
        freeData(&Hx_i);
        freeData(&Hy_i);
        freeData(&Ez_i);

        freeData(&Hx_s);
        freeData(&Hy_s);
        freeData(&Ez_s);
        freeData(&Ez_s_pre);

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

