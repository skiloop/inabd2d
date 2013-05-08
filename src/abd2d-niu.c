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
// 其他 --- 唉里
#define NIU_TYPE 5

//是否通过振幅算Eeff
//设为 1 算振幅
int if_erms_E_max=3; // if set niutype=4 then this is reset to 1

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
void input(int argc, char*argv[]);

int main(int argc, char *argv[]) {
    input(argc, argv);
    if(niutype==4){
        if_erms_E_max=1;
    }
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
#define MAX_BUFFER 201

void input(int argc, char*argv[]) {
    void PrintInput();
    void help();
    int i;
    //int cnt;
    //char buffer[MAX_BUFFER], *pstr;
    //double tdoub;
    if (argc < 2) {
        fprintf(stderr, "Invalid input!\n");
        help();
        exit(-1);
    }
    i = 1;
    while (i < argc) {
        if (strncmp(argv[i], "--niu-type=", 11) == 0) {
            niutype = atoi(argv[i] + 11);
            i++;
            continue;
        }
        if (strncmp(argv[i], "--matlab=1", 10) == 0) {
            IsMatlabSim = 1;
            i++;
            continue;
        }
        if (strncmp(argv[i], "--rei=", 6) == 0) {
            rei = atof(argv[i] + 6);
            i++;
            continue;
        }
        if (strncmp(argv[i], "--niu-e=", 8) == 0) {
            mu_e = atof(argv[i] + 8);
            i++;
            continue;
        }

        if (strncmp(argv[i], "--tm=1", 6) == 0) {
            IsTMx = 0;
            IsTEx = 1;
            i++;
            continue;
        }
        //MAXWELL_MESH_SIZE
        if (strncmp(argv[i], "--maxwell-grid=", 15) == 0) {
            m = atoi(argv[i] + 12);
            if (m < 4 || m > 40)m = MAXWELL_MESH_SIZE;
            i++;
            continue;
        }
        if (strncmp(argv[i], "--fine-grid=", 12) == 0) {
            m = atoi(argv[i] + 12);
            if (m < 4 || m > 40)m = FINE_GRID_SIZE;
            i++;
            continue;
        }
        if (strncmp(argv[i], "--total-time=", 13) == 0) {
            totaltime = atof(argv[i] + 13);
            if (totaltime <= 0)totaltime = TOTAL_TIME;
            i++;
            continue;
        }
        if (strncmp(argv[i], "--is-connect=", 13) == 0) {
            isConnect = atoi(argv[i] + 13);
            i++;
            continue;
        }
        if (strncmp(argv[i], "--with-density=", 15) == 0) {
            IfWithDensity = atoi(argv[i] + 15);
            i++;
            continue;
        }

        i++;
    }
    PrintInput();

}

void help() {
    printf("Usage:\n");
    printf("XXXX option\n\n");
    printf("--niu-type=n\tChoose niu formula type to compute.\n");
    printf("\t1\t------\tMorrow and Lowke's\n");
    printf("\t2\t------\tNikonov's\n");
    printf("\t3\t------\tKang's\n");
    printf("\tothers\t------\tdefault\n");
    printf("--matlab=[1,0]\tUse Matlab engine to show results\n");
    printf("--rei=XXX\tSet rei\n");
    printf("--niu-e=XXX\tSet Niu_e\n");
    printf("--tm=1\tset to TM mode\n");
    printf("--total-time=[>0]\ttotal simulation time\n");
    printf("--fine-grid=\thow many fine grid cells per Maxwell cell\n");
    printf("--maxwell-grid=\thow many Maxwell cells per wavelength\n");
    printf("--is-connect=[0,1]\tuse connecting interface or not\n");
    printf("--with-density=[0,1]\twether with density\n");
}

void PrintInput() {
    printf("--niu-type=%d\n", niutype);
    printf("--matlab=%d\n", IsMatlabSim);
    printf("--rei=%5.3e\n", rei);
    printf("--niu-e=%5.3e\n", mu_e);
    printf("--tm=%d\n", IsTMx);
    printf("--fine-grid=%d\n", m);
    printf("--total-time=%5.3e\n", totaltime);
    printf("--is-connect=%d\n", isConnect);
    printf("--with-density=%d\n", IfWithDensity);
}

