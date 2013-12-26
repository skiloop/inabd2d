
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "dataType.h"

MyDataF De;
MyDataF Da;
MyDataF D_kasi_max;
MyDataF mu_e;
MyDataF mu_i;
MyDataF vi, va;
MyDataF kasi;
MyDataF alpha;
MyDataF p = 760.0; //压强
MyDataF vm = 760 * 5.3e9;

MyDataF rei;

//FDTD DATA
MyDataF dt, dt_F, dt_M;
MyDataF half_dt;
MyDataF ds_F, ds_M;
MyDataF dx;
MyDataF dy;

int m = FINE_GRID_SIZE;
int MultiSize;
int n;

//DOMAIN DATA
int nx, nxp1, nxm1;
int ny, nyp1, nym1;
int nbound;

//PARAMETERS OF INCIDENT WAVE
MyDataF f; //frequency
MyDataF k; //
MyDataF T; //
MyDataF E_0 = E_MAX; // E field amptidute
MyDataF E0, H0;
MyDataF Hx0, Hz0, Hy0, Ez0, Ex0, Ey0;
MyDataF Ratio_x, Ratio_y;
MyDataF lamda;
MyDataF omega;
MyDataF phi = 0; //incidence wave inject angle on x-axis 

int IsTMx = _SOURCE_TMX_;
int IsTEx = _SOURCE_TEX_;

int TotalTimeStep;
int SaveTimeStep = SAVE_LEAP;
int Density_Time_Step = DEN_TIME_STEPS;
MyDataF CurTime = 0; //current time
int save_vi = 0;

// type of niu function
int niutype = NIU_TYPE;
//Matlab simulation
int IsMatlabSim = 0;
MyDataF totaltime = TOTAL_TIME; //in nane second
int isConnect = ISCONNECT;
int IfWithDensity = IF_WITH_DENSITY;
int maxwellGridSize = MAXWELL_MESH_SIZE;

//////////////////////////////////
// source waveform
//////////////////////////////////
int srcType = DEFAULT_SOURCE_TYPE;
MyDataF t_0; // time delay
MyDataF tauSquare; // tau^2 for Gaussian wave
MyDataF(*pSource)(MyDataF t);

//Total fields zone position
int tpis;
int tpjs;
int tpie;
int tpje;

int pis = 0;
int pie = 0;
int pjs = 0;
int pje = 0;

int if_erms_E_max = IF_ERMS_E_MAX; // if set niutype=4 then this is reset to 1

// number of threads for openmpi
#ifdef _OPENMP
int thread_count = DEFAULT_THREAD_COUNT;
#endif



