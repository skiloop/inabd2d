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

    /**
     * Initial connecting interface
     */
    void initconnect();

    /**
     * Apply connecting interface to E fields
     * @param t
     */
    void econnect(MyDataF t);

    /**
     * Apply connecting interface to M fields
     * @param t
     */
    void mconnect(MyDataF t);
    void end_connect();

    /**
     *     this function define where the total field domain is
     *     Note: pis,pjs,pie,pje must be defined and initialzed 
     *               before calling this function
     */
    void IntTtlFldDmnBnd();

    /**
     * Initial source
     */
    void initSource();

    /**
     * 
     * @param t
     * @return 
     */
    MyDataF Source(MyDataF t);

    /**
     * 
     * @param t
     * @return 
     */
    MyDataF GaussianSource(MyDataF t);

    /**
     * 
     * @param t
     * @return 
     */
    MyDataF SineSource(MyDataF t);

    /**
     * 
     * @param t
     * @return 
     */
    MyDataF CosineSource(MyDataF t);

#ifdef	__cplusplus
}
#endif

#endif	/* CONNECTINGITERFACE_H */

