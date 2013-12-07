/* 
 * File:   absortBound.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 上午11:19
 */

#ifndef PML_H
#define	PML_H

#ifdef	__cplusplus
extern "C" {
#endif


    void UpdMagFldForPML_TEz(MyStruct hx, MyStruct hy, MyStruct ez);
    void UpdMagFldForPML_TMz(MyStruct hz, MyStruct ex, MyStruct ey);
    void UpdEltFldForPML_TMz(MyStruct ex, MyStruct ey, MyStruct hz);
    void UpdEltFldForPML_TEz(MyStruct ez, MyStruct hx, MyStruct hy);
    void UpdateMFieldForPML(MyStruct hx, MyStruct hy, MyStruct hz, MyStruct ex, MyStruct ey, MyStruct ez);
    void UpdateEFieldForPML(MyStruct hx, MyStruct hy, MyStruct hz, MyStruct ex, MyStruct ey, MyStruct ez);
    void InitPMLBoundTEx();
    void InitPMLBoundTMx();
    void InitPMLBound();
    void FreePMLSpace();
    void InitBndCtrl();

#ifdef	__cplusplus
}
#endif

#endif	/* PML_H */

