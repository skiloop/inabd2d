/* 
 * File:   init.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 上午10:42
 */

#ifndef INITIAL_H
#define	INITIAL_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "dataType.h"
    
    // field create and zone difenition function
    void CalDomainSize() ;
    
    // field initial functions
void initCommonData();
void Init_ne();
void InitBeta(MyStruct* stru, MyDataF gamma);
void InitCeh(MyStruct beta);
void InitCev(MyStruct beta, MyDataF alpha_ev);
void InitVelocityCoeff_Equatioan_Five() ;

// funtion with no density
void InitCeeNoDen();
void InitCehNoDen();
void InitCevNoDen() ;

// coefficient initial and update function
void InitCoeff();
void UpdateCoeff();

//
void FreeSpace();

#ifdef	__cplusplus
}
#endif

#endif	/* INITIAL_H */

