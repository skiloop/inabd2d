/* 
 * File:   fdtd.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 上午11:08
 */

#ifndef FDTD_H
#define	FDTD_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "dataType.h"

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
    //Coeffients    
    extern MyDataF Chxez, Chyez, Chzex, Chzey;
    extern MyStruct Cevx, Cevy, Cevz;
    extern MyStruct Ceex, Ceey, Ceez;
    extern MyStruct Cehx, Cehy, Cehz;
    extern MyStruct Cvvx, Cvvy, Cvvz;
    extern MyStruct Cvex, Cvey, Cvez;
    extern MyDataF Cve;

    void UpdateEField();
    void UpdateMField();
    void fdtd();

    void my_pause();
    void UpdateDensity();
    void InterpolatErms();
    void DensityBound(MyStruct stru, int bndwidth, const int swidth);
    void CalErmsAtCoarseGrid();
    void CalSumESqrt();
    void CalSumESqrt_Emax();
    void calErmsAtCoarsePoint();
    void calErmsAtCoarsePoint_Max();

#ifdef	__cplusplus
}
#endif

#endif	/* FDTD_H */

