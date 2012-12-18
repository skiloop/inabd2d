#define _CRT_SECURE_NO_DEPRECATE

#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<string.h>
#include<time.h>

#ifdef _MSC_VER
#include <float.h>
#define isnan _isnan
#endif

/////////////////////////////////
//FDTD DEFINITION
/////////////////////////////////
#define MAX_LIMIT 1e80
#define NE_MAX_LIMIT 1E307
#define COURANT_FACTOR 0.5
#define CFL_FACTOR 0.4
#define MAXWELL_MESH_SIZE 1.0/50.0
#define NUMBER_OF_CELLS_IN_PML_BOUND 20
#define SCATTER_FIELD_DOMAIN_BND_SIZE 5
#define NUMBER_OF_WAVELENGTHS_IN_DOMAIN 3
#define DEN_TIME_STEPS 2000000

#define TOTAL_TIME 150e-9 // in second
#define SAVE_LEAP 150	//sample electric field components and ne every SAVE_LEAP steps
#define SAVE_ERMS_LEAP 110
#define PULSE_LENGGTH_IN_TIME
#define FINE_GRID_SIZE 10
#define SCATTERED_FORMATION
#define MATLAB_SIMULATION_
#define _SOURCE_TEX_ 0
#define _SOURCE_TMX_ 1 

#define E_0 6e6
#define NE0 1e13

#define FREQUENCY 110E9
#define MAX_NE 1E29
#define INC_ANGLE 0.0*M_PI
//////////////////////////////////////
//display definition
/////////////////////////////////////
#define FIELD_TO_DISPLAY0 Ex_s
#define FIELD_TO_DISPLAY1 Ey_s
#define FIELD_TO_DISPLAY2_ Hz_s
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



#include"data_type_definetion.h"
#include"common_data.h"
#include"init_com_data.h"
#include"bndctrl.h"
#include"init_data.h"
#include"calculate_domain_size.h"


#include"test.h"
#include"boundary.h"
#include"source.h"
#include"scatter_total_bnd.h"
#include"updatecoef.h"
/*******************************************/


#include"simulate.h"
/*******************************************/

#include"scatfdtd.h"
#include"savedata.h"
#include"freespace.h"

int main(int argc,char *argv[]) {
    InitComData();
    CalDomainSize();
    InitBndCtrl();
    InitCoeff();
    InitPMLBound();
    //InitPMLBoundi();
    IntTtlFldDmnBnd();
    //InitIncFdtd();
    scatfdtd();
    FreeDelayArrays();
    FreePMLSpace();
    SaveToFile();
    FreeSpace();
    //FreePMLSpace();
    return 0;
}
