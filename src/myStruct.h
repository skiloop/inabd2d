/* 
 * File:   myStruct.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 上午10:57
 */

#ifndef MYSTRUCT_H
#define	MYSTRUCT_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct {
        int nx;
        int ny;
        MyDataF* data;
    } MyStruct;

    /**
     * 
     * @param mst
     */
    void ResetStructData(MyStruct mst);
    /**
     * 
     * @param mst
     * @param val
     */
    void ResetStructDataR(MyStruct mst, MyDataF val);
    /**
     * 
     * @param nnx
     * @param nny
     * @param ms
     */
    void InitMyStr(int nnx, int nny, MyStruct* ms);
    /**
     * 
     * @param mst
     */
    void ResetStructData2(MyStruct mst);
    /**
     * 
     * @param mst
     */
    void PrintData(MyStruct mst);
    /**
     * 
     * @param stru
     */
    void CheckStruct(const MyStruct stru);
    /**
     * 
     * @param num
     * @param fname
     * @param data
     */
    void CaptData(const int num, const char * const fname, const MyStruct data);
    /**
     * 
     * @param num
     * @param fname
     * @param data
     * @param is
     * @param ie
     * @param js
     * @param je
     */
    void CaptDataNoPML(const int num, const char * const fname, const MyStruct data, int is, int ie, int js, int je);
    /**
     * 
     * @param num
     * @param fname
     * @param data
     * @param jcenter
     * @param is
     * @param ie
     */
    void CaptDataCenter(const int num, const char * const fname, const MyStruct data, int jcenter, int is, int ie);
    /**
     * 
     * @param num
     * @param fname
     * @param data
     * @param p
     */
    void CaptDataM(const int num, const char*fname, const MyStruct data, int p);
    /**
     * 
     * @param num
     * @param fname
     * @param data
     * @param p
     * @param is
     * @param ie
     * @param js
     * @param je
     */
    void CaptDataMNoPML(const int num, const char*fname, const MyStruct data, int p, int is, int ie, int js, int je);
    /**
     * 
     */
    void SaveToFile(); //defined in savedata.h
    /**
     * 
     * @param coo
     * @param nx
     * @param ny
     * @param InitValue
     */
    void initMyStruct(MyStruct *coo, int nx, int ny, MyDataF InitValue);

    /**
     * 
     * @param st
     * @param stpre
     * @return 
     */
    int BackupMyStruct(MyStruct st, MyStruct stpre) ;
#ifdef	__cplusplus
}
#endif

#endif	/* MYSTRUCT_H */


