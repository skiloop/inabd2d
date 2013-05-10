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
    static const MyDataF mu_0 = 1.257e-6;
    static const MyDataF eps_0 = 8.854e-12;
    static const MyDataF length = 1.0;
    static const MyDataF e = 1.602e-19; //electron charge
    static const MyDataF c = 2.998E8; //speed of light
    static const MyDataF me = 9.110e-31; //electron mass
    static const MyDataF CFL_factor = CFL_FACTOR;
    static const MyDataF CourantFactor = COURANT_FACTOR;

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
    extern MyDataF E0, H0;
    extern MyDataF Hx0, Hz0, Hy0, Ez0, Ex0, Ey0;
    extern MyDataF Ratio_x, Ratio_y;
    extern MyDataF lamda;
    extern MyDataF omega;
    extern MyDataF phi; //incidence wave inject angle on x-axis 


    //EM field
    extern MyStruct Ey_i, Ex_i, Ez_i, Hx_i, Hy_i, Hz_i;
    extern MyStruct Ey_s, Ex_s, Ez_s, Hx_s, Hy_s, Hz_s;
    extern MyStruct Ey_s_pre, Ex_s_pre, Ez_s_pre;
    extern MyStruct Vey, Vex, Vez;
    extern MyStruct ne, ne_pre;
    extern MyStruct beta;
    extern MyStruct Erms;
    extern MyStruct Ermsx, Ermsy;
    extern MyStruct Nu_c; // 循环利用的碰撞率
    //MyStruct Nu_c_pre;//上一步的碰撞率（循环利用的碰撞率）
    //Coeffients
    /*******************************/
    extern MyDataF Chxez, Chyez, Chzex, Chzey;

    extern MyStruct Cevx, Cevy, Cevz;
    extern MyStruct Ceex, Ceey, Ceez;
    extern MyStruct Cehx, Cehy, Cehz;

    extern MyDataF Cve;

    extern MyStruct Cvvx, Cvvy, Cvvz;
    extern MyStruct Cvex, Cvey, Cvez;

    /*******************************/

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
    extern FILE *filedat;
    MyDataF totaltime; //in nane second
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

#ifdef	__cplusplus
}
#endif

#endif	/* COMMONDATA_H */

