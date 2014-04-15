/* 
 * File:   commonData.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 上午10:24
 */

#ifndef COMMONDATA_H
#define	COMMONDATA_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "dataType.h"

    // common constants    
    static const MyDataF mu_0 = 1.2566370614359173e-06;
    static const MyDataF eps_0 = 8.8541878176203892e-12;
    static const MyDataF length = 1.0;
    static const MyDataF e = 1.602e-19; //electron charge
    static const MyDataF c = 299792458.0; //speed of light
    static const MyDataF me = 9.110e-31; //electron mass
    static const MyDataF CFL_factor = CFL_FACTOR;
    static const MyDataF CourantFactor = COURANT_FACTOR;
    static const MyDataF NeutralGasDensityCM = 2.44e19; // neutral gas number density in cm^-3
    static const MyDataF NeutralGasDensityM = 2.44e25; // neutral gas number density in m^-3

    extern MyDataF De;
    extern MyDataF Da;
    extern MyDataF D_kasi_max;
    extern MyDataF mu_e;
    extern MyDataF mu_i;
    extern MyDataF vi, va;
    extern MyDataF kasi;
    extern MyDataF alpha;
    extern MyDataF p; //压强

    extern MyDataF vm;

    extern MyDataF rei;

    //FDTD DATA
    extern MyDataF dt, dt_F, dt_M;
    extern MyDataF half_dt;
    extern MyDataF ds_F, ds_M;
    extern MyDataF dx;
    extern MyDataF dy;
    extern int m;
    extern int n;
    extern int MultiSize;

    //DOMAIN DATA
    extern int nx, nxp1, nxm1;
    extern int ny, nyp1, nym1;
    extern int nbound;

    //PARAMETERS OF INCIDENT WAVE
    extern MyDataF f; //frequency
    extern MyDataF k; //
    extern MyDataF T; //
    extern MyDataF E0, H0, E_0;
    extern MyDataF RatioHx, RatioHz, RatioHy, RatioEz, RatioEx, RatioEy;

    extern MyDataF lamda;
    extern MyDataF omega;
    extern MyDataF phi; //incidence wave inject angle on x-axis 
    extern MyDataF cos_phi, sin_phi;
    extern MyDataF psi; //phase angel
    extern MyDataF cos_psi, sin_psi;

    extern int IsTMx;
    extern int IsTEx;

    extern int TotalTimeStep;
    extern int SaveTimeStep;
    extern int Density_Time_Step;
    extern MyDataF CurTime; //current time
    extern int save_vi;

    // type of niu function
    extern int niutype;
    //Matlab simulation
    extern int IsMatlabSim;
    extern MyDataF totaltime; //in nane second
    extern int isConnect;
    extern int IfWithDensity;
    extern int maxwellGridSize;

    //Total fields zone position
    extern int tpis;
    extern int tpjs;
    extern int tpie;
    extern int tpje;

    //absorb bound fields zone position
    extern int pis;
    extern int pjs;
    extern int pie;
    extern int pje;

    extern int if_erms_E_max; // if set niutype=4 then this is reset to 1
#ifdef _OPENMP
    extern int thread_count;
#endif

    // source type
    extern int srcType;
    extern MyDataF t_0; // time delay
    extern MyDataF tauSquare; // tau^2 for Gaussian wave
    extern MyDataF(*pSource)(MyDataF t);

    // em fields
    extern MyStruct Ex, Ey, Ez, Hx, Hy, Hz, Ex_pre, Ey_pre, Ez_pre;

#ifdef	__cplusplus
}
#endif

#endif	/* COMMONDATA_H */

