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

void UpdateEField() ;
void UpdateMField() ;
void fdtd();

void my_pause() ;
void UpdateDensity() ;
void InterpolatErms() ;
void DensityBound(MyStruct stru, int bndwidth, const int swidth) ;
void CalErmsAtCoarseGrid() ;
void CalSumESqrt() ;
void CalSumESqrt_Emax() ;
void calErmsAtCoarsePoint() ;
void calErmsAtCoarsePoint_Max() ;

#ifdef	__cplusplus
}
#endif

#endif	/* FDTD_H */

