/* 
 * File:   common.h
 * Author: skiloop
 *
 * Created on 2013年5月9日, 上午10:20
 */

#ifndef COMMON_H
#define	COMMON_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <math.h>

#ifdef _MSC_VER
#include <float.h>
#define isnan _isnan
#endif

#ifndef M_PI
#define M_PI 3.1415926
#endif
    /////////////////////////////////
    //FDTD DEFINITION
    /////////////////////////////////
#define MAX_LIMIT 1e80
#define NE_MAX_LIMIT 1E307
#define COURANT_FACTOR 0.5
#define CFL_FACTOR 0.4
    // maxwll 网格大小
#define MAXWELL_MESH_SIZE 50
    //PML网格数
#define NUMBER_OF_CELLS_IN_PML_BOUND 20
    //
#define SCATTER_FIELD_DOMAIN_BND_SIZE 5
#define NUMBER_OF_WAVELENGTHS_IN_DOMAIN 3
#define DEN_TIME_STEPS 2000000

#define TOTAL_TIME 100.0 // in nane second
    //  间隔多少网格保存
#define SAVE_LEAP 50	//sample electric field components and ne every SAVE_LEAP steps
    // 
#define SAVE_ERMS_LEAP 110
#define PULSE_LENGGTH_IN_TIME
    // 亚网格相对大小：一个Maxwell有多少个亚网格
#define FINE_GRID_SIZE 10
#define SCATTERED_FORMATION

    // 波的模式 
    // TEx --- Ez Hx Hy
    // TMx --- Hz Ex Ey
    // 单模式的话要设置下面两个一个为0另一个为1
    // 如选择 TEx 设置 _SOURCE_TEX_ 1
    // _SOURCE_TMX_ 0
#define _SOURCE_TEX_ 1
#define _SOURCE_TMX_ 0 

    // 电离公式
    // 1 ---- MorrowAndLowke
    // 2 ---- Nikonov
    // 3 ---- Kang
    // 4 ---- Erms as Emax with Zhao
    // 5 ---- Zhao
    // 其他 --- Ali
#define NIU_TYPE 5

    //是否通过振幅算Eeff
    //设为 1 算振幅
   //const  int if_erms_E_max = 3; // if set niutype=4 then this is reset to 1
    // 在abd2d.c里改
#define IF_ERMS_E_MAX 3 // default values for if_erms_E_max

#define E_0 6.0e6
#define NE0 1e13

#define FREQUENCY 110E9
#define MAX_NE 1E24
#define INC_ANGLE 0.0*M_PI
    //////////////////////////////////////
    //display definition
    /////////////////////////////////////
#define FIELD_TO_DISPLAY0 Ex_s
#define FIELD_TO_DISPLAY1 Ey_s
#define FIELD_TO_DISPLAY2 Hz_s

#define FIELD_TO_DISPLAY3 Hx_s
#define FIELD_TO_DISPLAY4 Hy_s
#define FIELD_TO_DISPLAY5 Ez_s

#define DISPLAY_NE


    ////////////////////////////////////
    //SAMPLE DEFINITION
    ///////////////////////////////////
#define LEAPSTEP_OF_DISPLAY 1
#define DTF_VI_LIMIT 2.0e10
#define LEAPSTEP_OF_CAPTURE 40
#define CAPTURE_FIELD ne

    ///////////////////////////
    //SOURCES POSITION
    ///////////////////////////
#define CELLS_INSIDE_I 0
#define CELLS_INSIDE_S 0
#define SOURCE_POS_IN_Y 0
#define SOURCE_POS_IN_X 0
#define SOURCES_SIZE 1
#define ISCONNECT 1
#define IF_WITH_DENSITY 1

	////////////////////////////
	// OUTPUT FILE TYPE
	////////////////////////////
#define ORIGN_TYPE (1) // output for orgin
#define MATLAB_TYPE (2) // output for matlab
#define SAVE_TYPE (MATLAB_TYPE)


#ifdef	__cplusplus
}
#endif

#endif	/* COMMON_H */

