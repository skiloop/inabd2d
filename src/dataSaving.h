/* 
 * File:   dataSaving.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 下午4:25
 */

#ifndef DATASAVING_H
#define	DATASAVING_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "dataType.h"

    extern MyDataF *srcdat;
    void CheckIfOverLimit(int i, int j, MyDataF data);
    void saveEfieldCenterTMz(int timestep);
    void CheckIfNeOverLimit(int i, int j, MyDataF data);
    void chkgrad(MyDataF str, MyDataF strpre, int i, int j);
    void CheckIfNonZeros(MyStruct stru);
    void SaveCapField(const int timestep);
    void SaveErms(int nestep);
    void SaveRows(const MyStruct stru, int row, char *fname);
    void calsampos(int *xpos, int *ypos);
    void InitCapture(int ttstep);
    void CaptureE(MyDataF *CapEf, int CurTimeStep, int i, int j);
    void scpfld(MyDataF *CapEf, int ttstep);
    void SaveData(MyStruct data, char* filename, FILE *hfile);
    void SaveToFile();


#ifdef	__cplusplus
}
#endif

#endif	/* DATASAVING_H */

