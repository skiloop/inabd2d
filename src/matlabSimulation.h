/* 
 * File:   matlabSimulation.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 上午11:19
 */

#ifndef MATLABSIMULATION_H
#define	MATLABSIMULATION_H

#ifdef	__cplusplus
extern "C" {
#endif

    //    
    //#ifdef MATLAB_SIMULATION
    //    
    ////#include"matrix.h"
    //#include"engine.h"
    //#include"mex.h"
    //    
    //void PreSimuAdd1(Engine *ep, const MyStruct stru) ;
    //void PreSimuAdd0(Engine *ep, const MyStruct stru) ;
    //void PreSimuAdd2(Engine *ep, const MyStruct stru) ;
    //void Simulate0(Engine *ep, const MyStruct stru) ;
    //void Simulate1(Engine *ep, const MyStruct stru) ;
    //void Simulate2(Engine *ep, const MyStruct stru) ;
    //void PreSimuAddNE(Engine *ep, const MyStruct stru) ;
    //void SimulateNe(Engine *ep, const MyStruct stru) ;
    //
    //#endif
    void DispEMFields(const int timestep);
    void EndMatlabSim();
    void DispNe(const int timestep);
    void InitSim();
#ifdef	__cplusplus
}
#endif

#endif	/* MATLABSIMULATION_H */


