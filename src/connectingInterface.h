/* 
 * File:   connectingIterface.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 上午11:18
 */

#ifndef CONNECTINGITERFACE_H
#define	CONNECTINGITERFACE_H

#ifdef	__cplusplus
extern "C" {
#endif

    void initconnect();
    void econnect(MyDataF t);
    void mconnect(MyDataF t);
    void end_connect();

    void IntTtlFldDmnBnd();
    void Adjust_E_Field(MyDataF pre_t); //adjust TEx electric field which at total fields boundary
    void Adjust_M_Field(MyDataF pre_t);
    void FreeDelayArrays();

    MyDataF Source(MyDataF t);


#ifdef	__cplusplus
}
#endif

#endif	/* CONNECTINGITERFACE_H */

