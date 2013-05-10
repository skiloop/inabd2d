/* 
 * File:   fluid.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 上午11:11
 */

#ifndef FLUID_H
#define	FLUID_H

#ifdef	__cplusplus
extern "C" {
#endif

//velocity functions
void UpdateVelocity();

// density functions
void UpdateDensity() ;
void DensityBound(MyStruct stru, int bndwidth, const int swidth) ;

// Erms and Eeff functions
void CalErmsAtCoarseGrid() ;
void CalSumESqrt() ;
void CalSumESqrt_Emax() ;
void InterpolatErms() ;
void calErmsAtCoarsePoint() ;
void calErmsAtCoarsePoint_Max() ;

void my_pause() ;
#ifdef	__cplusplus
}
#endif

#endif	/* FLUID_H */

