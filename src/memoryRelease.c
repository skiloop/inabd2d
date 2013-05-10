
#include "common.h"
#include "dataType.h"
#include "fdtd.h"
#include "commonData.h"

void FreeSpace() {
    if (IsTMx) {
        freeData(&Ex_s);
        freeData(&Ex_s_pre);

        freeData(&Vex);
        freeData(&Cevx);
        freeData(&Ceex);
        freeData(&Cehx);

        freeData(&Ey_s);
        freeData(&Ey_s_pre);

        freeData(&Vey);
        freeData(&Cevy);
        freeData(&Ceey);
        freeData(&Cehy);
        freeData(&Hz_s);
    }
    if (IsTEx) {
        freeData(&Hx_i);
        freeData(&Hy_i);
        freeData(&Ez_i);

        freeData(&Hx_s);
        freeData(&Hy_s);
        freeData(&Ez_s);
        freeData(&Ez_s_pre);

        freeData(&Vez);
        freeData(&Ceez);
        freeData(&Cevz);
        freeData(&Cehz);
    }
    freeData(&Cvvx);
    freeData(&Cvvy);
    freeData(&Cvvz);
    freeData(&Cvex);
    freeData(&Cvey);
    freeData(&Cvez);
    freeData(&ne);
    freeData(&Erms);
    freeData(&ne_pre);
    freeData(&beta);
    freeData(&Nu_c);
}

