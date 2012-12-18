#ifndef COMMOM_DATA_INCLUDE
#define COMMOM_DATA_INCLUDE
#pragma once
#include"data_type_definetion.h"


MyDataF De;
MyDataF Da;
MyDataF D_kasi_max;
MyDataF mu_e;
MyDataF mu_i;
MyDataF vi,va;
MyDataF kasi;
MyDataF alpha;
MyDataF p	=	760.0;
CMyDataF mu_0	=	1.257e-6;
CMyDataF eps_0	=	8.854e-12;
CMyDataF vm	=	760*5.3e9;
CMyDataF length	=	1.0;


//electron charge
CMyDataF e=1.602e-19;

//speed of light
CMyDataF c=2.998E8;

//electron mass
CMyDataF me=9.110e-31;
MyDataF rei;


//FDTD DATA
MyDataF dt,dt_F,dt_M;
MyDataF half_dt;
MyDataF ds_F,ds_M;
MyDataF dx;
MyDataF dy;
CMyDataF CFL_factor = CFL_FACTOR;
CMyDataF CourantFactor = COURANT_FACTOR;
int m = FINE_GRID_SIZE;
int n;

//DOMAIN DATA
int nx,nxp1,nxm1;
int ny,nyp1,nym1;
int nbound;

//PARAMETERS OF INCIDENT WAVE
MyDataF f;//frequency
MyDataF k;//
MyDataF T;//
MyDataF E0,H0;
MyDataF Hx0,Hz0,Hy0,Ez0,Ex0,Ey0;
MyDataF Ratio_x,Ratio_y;
MyDataF lamda;
MyDataF omega;
MyDataF phi=INC_ANGLE;//incidence wave inject angle on x-axis 


//EM field
MyStruct Ey_i,Ex_i,Ez_i,Hx_i,Hy_i,Hz_i;
MyStruct Ey_s,Ex_s,Ez_s,Hx_s,Hy_s,Hz_s;
MyStruct Ey_s_pre,Ex_s_pre,Ez_s_pre;
MyStruct Vey,Vex,Vez;
MyStruct ne,ne_pre;
MyStruct beta;
MyStruct Erms;
MyStruct Ermsx,Ermsy;
//Coeffients
/*******************************/
MyDataF Chxez,Chyez,Chzex,Chzey;

MyStruct Cevx,Cevy,Cevz;
MyStruct Ceex,Ceey,Ceez;
MyStruct Cehx,Cehy,Cehz;

MyDataF Cve;
/*******************************/

int IsTMx = _SOURCE_TMX_;
int IsTEx = _SOURCE_TEX_;

int TotalTimeStep;
int SaveTimeStep=SAVE_LEAP;
int Density_Time_Step = DEN_TIME_STEPS;
MyDataF CurTime=0;//current time
int save_vi = 0;

// type of niu function
int niutype = 4;
//Matlab simulation
int IsMatlabSim = 0;
FILE *filedat;
double totaltime = TOTAL_TIME;//in nane second
int isConnect = ISCONNECT;
int IfWithDensity = IF_WITH_DENSITY;
int maxwellGridSize=MAXWELL_MESH_SIZE;
#endif

