
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dataType.h"
#include "commonData.h"

MyStruct Hzx_xn;
MyStruct Hzx_xp;

MyStruct Hzx_yn;
MyStruct Hzx_yp;

MyStruct Hzy_xn;
MyStruct Hzy_xp;

MyStruct Hzy_yn;
MyStruct Hzy_yp;


MyStruct Ceye_xn;
MyStruct Ceyhz_xn;
MyStruct Chzxh_xn;
MyStruct Chzyh_xn;
MyStruct Chzyex_xn;
MyStruct Chzxey_xn;

MyStruct Ceye_xp;
MyStruct Ceyhz_xp;
MyStruct Chzxh_xp;
MyStruct Chzyh_xp;
MyStruct Chzyex_xp;
MyStruct Chzxey_xp;

MyStruct Cexe_yp;
MyStruct Cexhz_yp;
MyStruct Chzxh_yp;
MyStruct Chzyh_yp;
MyStruct Chzxey_yp;
MyStruct Chzyex_yp;

MyStruct Cexe_yn;
MyStruct Cexhz_yn;
MyStruct Chzxh_yn;
MyStruct Chzyh_yn;
MyStruct Chzxey_yn;
MyStruct Chzyex_yn;




MyStruct Ezx_xn;
MyStruct Ezx_xp;

MyStruct Ezx_yn;
MyStruct Ezx_yp;

MyStruct Ezy_xn;
MyStruct Ezy_xp;

MyStruct Ezy_yn;
MyStruct Ezy_yp;

MyStruct Cezxe_xn;
MyStruct Cezxhy_xn;
MyStruct Chyh_xn;
MyStruct Chyez_xn;
MyStruct Cezye_xn;
MyStruct Cezyhx_xn;

MyStruct Cezye_yn;
MyStruct Cezyhx_yn;
MyStruct Cezxe_yn;
MyStruct Cezxhy_yn;
MyStruct Chxh_yn;
MyStruct Chxez_yn;


MyStruct Cezxe_xp;
MyStruct Cezxhy_xp;
MyStruct Cezye_xp;
MyStruct Cezyhx_xp;
MyStruct Chyh_xp;
MyStruct Chyez_xp;

MyStruct Cezxe_yp;
MyStruct Cezxhy_yp;
MyStruct Cezye_yp;
MyStruct Cezyhx_yp;
MyStruct Chxh_yp;
MyStruct Chxez_yp;

int IsPMLxn = 0;
int IsPMLxp = 0;
int IsPMLyn = 0;
int IsPMLyp = 0;

int IsAnySidePML = 0;

int PMLCellxn = 0;
int PMLCellxp = 0;
int PMLCellyn = 0;
int PMLCellyp = 0;

int pml_order = 2;
MyDataF pml_R_0 = 1e-8;


void InitBndCtrl() {
    IsAnySidePML = 1;
    PMLCellxp = nbound;
    PMLCellyp = nbound;
    PMLCellxn = nbound;
    PMLCellyn = nbound;
    IsPMLxp = 1;
    IsPMLyp = 1;
    IsPMLxn = 1;
    IsPMLyn = 1;

    pis = PMLCellxn;
    pie = nx - PMLCellxp;

    pjs = PMLCellyn;
    pje = ny - PMLCellyp;
}

/*******************************************************************************************************
UpdMagFldForPML_TEz.h
update magnetic fields at PML regions
TEz
 ***************************************************************************************************/
void UpdMagFldForPML_TEz(MyStruct hx, MyStruct hy, MyStruct ez) {
    int i, j;
    if (IsPMLxn) {
        for (i = 0; i < pis; i++)
            for (j = 1; j <= nym1; j++) {
                hy.data[i * hy.ny + j] = Chyh_xn.data[i * Chyh_xn.ny + j - 1] * hy.data[i * hy.ny + j]
                        + Chyez_xn.data[i * Chyez_xn.ny + j - 1]*(ez.data[(i + 1) * ez.ny + j] - ez.data[i * ez.ny + j]);
                /*******************************************************************************************/
            }
    }
    if (IsPMLxp) {
        for (i = pie; i <= nxm1; i++)
            for (j = 1; j <= nym1; j++) {
                hy.data[i * hy.ny + j] = Chyh_xp.data[(i - pie) * Chyh_xp.ny + j - 1] * hy.data[i * hy.ny + j]
                        + Chyez_xp.data[(i - pie) * Chyez_xp.ny + j - 1]*(ez.data[(i + 1) * ez.ny + j] - ez.data[i * ez.ny + j]);
                /*******************************************************************************************/
            }
    }

    if (IsPMLyn) {
        for (i = 1; i <= nxm1; i++)
            for (j = 0; j < pjs; j++)
                hx.data[i * hx.ny + j] = Chxh_yn.data[(i - 1) * Chxh_yn.ny + j] * hx.data[i * hx.ny + j]
                    + Chxez_yn.data[(i - 1) * Chxez_yn.ny + j]*(ez.data[i * ez.ny + j + 1] - ez.data[i * ez.ny + j]);
    }

    if (IsPMLyp) {
        for (i = 1; i <= nxm1; i++)
            for (j = pje; j <= nym1; j++)
                hx.data[i * hx.ny + j] = Chxh_yp.data[(i - 1) * Chxh_yp.ny + j - pje] * hx.data[i * hx.ny + j]
                    + Chxez_yp.data[(i - 1) * Chxez_yp.ny + j - pje]*(ez.data[i * ez.ny + j + 1] - ez.data[i * ez.ny + j]);
    }
}

/*********************************************************************
UpdMagFldForPML_TMz.h
update magnetic fields at PML regions
TMz
 *********************************************************************/
void UpdMagFldForPML_TMz(MyStruct hz, MyStruct ex, MyStruct ey) {
    int i, j;
    if (IsPMLxn) {
        for (i = 0; i < Hzx_xn.nx; ++i)
            for (j = 0; j < Hzx_xn.ny; j++)
                Hzx_xn.data[i * Hzx_xn.ny + j] = Chzxh_xn.data[i * Chzxh_xn.ny + j] * Hzx_xn.data[i * Hzx_xn.ny + j]
                    + Chzxey_xn.data[i * Chzxey_xn.ny + j]*(ey.data[(i + 1) * ey.ny + j] - ey.data[i * ey.ny + j]);
        for (i = 0; i < Hzy_xn.nx; ++i)
            for (j = 0; j < Hzy_xn.ny; j++) {
                Hzy_xn.data[i * Hzy_xn.ny + j] = Chzyh_xn.data[i * Chzyh_xn.ny + j] * Hzy_xn.data[i * Hzy_xn.ny + j]
                        + Chzyex_xn.data[i * Chzyex_xn.ny + j]*(ex.data[i * ex.ny + j + pjs + 1] - ex.data[i * ex.ny + j + pjs]);

            }
    }
    /*================================================================================================*/
    if (IsPMLxp) {
        for (i = 0; i < Hzx_xp.nx; ++i)
            for (j = 0; j < Hzx_xp.ny; j++) {
                Hzx_xp.data[i * Hzx_xp.ny + j] = Chzxh_xp.data[i * Chzxh_xp.ny + j] * Hzx_xp.data[i * Hzx_xp.ny + j]
                        + Chzxey_xp.data[i * Chzxey_xp.ny + j]*(ey.data[(i + pie + 1) * ey.ny + j] - ey.data[(i + pie) * ey.ny + j]);

                //if(fabs(Hzx_xp.data[i*Hzx_xp.ny+j])>(1.8e-2)*eps)
                //	puts("Hello!");
            }
        for (i = 0; i < Hzy_xp.nx; ++i)
            for (j = 0; j < Hzy_xp.ny; j++) {
                Hzy_xp.data[i * Hzy_xp.ny + j] = Chzyh_xp.data[i * Chzyh_xp.ny + j] * Hzy_xp.data[i * Hzy_xp.ny + j]
                        + Chzyex_xp.data[i * Chzyex_xp.ny + j]*(ex.data[(pie + i) * ex.ny + j + pjs + 1] - ex.data[(pie + i) * ex.ny + j + pjs]);
                //if(fabs(Hzx_xp.data[i*Hzx_xp.ny+j])>(1.8e-2)*eps)
                //	puts("Hello!");
            }
    }
    /*================================================================================================*/
    if (IsPMLyn) {
        for (i = 0; i < Hzx_yn.nx; ++i)
            for (j = 0; j < Hzx_yn.ny; j++)
                Hzx_yn.data[i * Hzx_yn.ny + j] = Chzxh_yn.data[i * Chzxh_yn.ny + j] * Hzx_yn.data[i * Hzx_yn.ny + j]
                    + Chzxey_yn.data[i * Chzxey_yn.ny + j]*(ey.data[(i + 1 + pis) * ey.ny + j] - ey.data[(i + pis) * ey.ny + j]);
        for (i = 0; i < Hzy_yn.nx; ++i)
            for (j = 0; j < Hzy_yn.ny; j++)
                Hzy_yn.data[i * Hzy_yn.ny + j] = Chzyh_yn.data[i * Chzyh_yn.ny + j] * Hzy_yn.data[i * Hzy_yn.ny + j]
                    + Chzyex_yn.data[i * Chzyex_yn.ny + j]*(ex.data[i * ex.ny + j + 1] - ex.data[i * ex.ny + j]);

    }
    /*===============================================================================================*/
    if (IsPMLyp) {
        for (i = 0; i < Hzx_yp.nx; ++i)
            for (j = 0; j < Hzx_yp.ny; j++) {
                Hzx_yp.data[i * Hzx_yp.ny + j] = Chzxh_yp.data[i * Chzxh_yp.ny + j] * Hzx_yp.data[i * Hzx_yp.ny + j]
                        + Chzxey_yp.data[i * Chzxey_yp.ny + j]*(ey.data[(i + pis + 1) * ey.ny + j + pje] - ey.data[(i + pis) * ey.ny + j + pje]);
                //if(fabs(Hzx_yp.data[i*Hzx_yp.ny+j])>((1e-2)*eps))
                //	i=i;
            }
        for (i = 0; i < Hzy_yp.nx; ++i)
            for (j = 0; j < Hzy_yp.ny; j++) {
                Hzy_yp.data[i * Hzy_yp.ny + j] = Chzyh_yp.data[i * Chzyh_yp.ny + j] * Hzy_yp.data[i * Hzy_yp.ny + j]
                        + Chzyex_yp.data[i * Chzyex_yp.ny + j]*(ex.data[i * ex.ny + j + pje + 1] - ex.data[i * ex.ny + j + pje]);
                //if(fabs(Hzy_yp.data[i*Hzy_yp.ny+j])>((1e-2)*eps))
                //	i=i;
            }
    }
    /*===============================================================================================*/
    //left bottom
    for (i = 0; i < pis; i++) {
        for (j = 0; j < pjs; j++) {
            hz.data[i * hz.ny + j] = Hzx_xn.data[i * Hzx_xn.ny + j] + Hzy_yn.data[i * Hzy_yn.ny + j];
        }
    }
    //left top
    for (i = 0; i < pis; i++) {
        for (j = pje; j < ny; j++) {
            hz.data[i * hz.ny + j] = Hzx_xn.data[i * Hzx_xn.ny + j] + Hzy_yp.data[i * Hzy_yp.ny + j - pje];
        }
    }
    //right bottom 
    for (i = pie; i < nx; i++) {
        for (j = 0; j < pjs; j++) {
            hz.data[i * hz.ny + j] = Hzx_xp.data[(i - pie) * Hzx_xp.ny + j] + Hzy_yn.data[i * Hzy_yn.ny + j];
        }
    }
    //right top
    for (i = pie; i < nx; i++) {
        for (j = pje; j < ny; j++) {
            hz.data[i * hz.ny + j] = Hzx_xp.data[(i - pie) * Hzx_xp.ny + j] + Hzy_yp.data[i * Hzy_yp.ny + j - pje];
        }
    }
    /*===============================================================================================*/
    for (i = 0; i < pis; i++) {
        for (j = pjs; j < pje; j++) {
            hz.data[i * hz.ny + j] = Hzx_xn.data[i * Hzx_xn.ny + j] + Hzy_xn.data[i * Hzy_xn.ny + j - pjs];
        }
    }
    for (i = 0; i < PMLCellxp; i++) {
        for (j = pjs; j < pje; j++) {
            hz.data[(i + pie) * hz.ny + j] = Hzx_xp.data[i * Hzx_xp.ny + j] + Hzy_xp.data[i * Hzy_xp.ny + j - pjs];
        }
    }
    for (i = pis; i < pie; i++) {
        for (j = 0; j < pjs; j++) {
            hz.data[i * hz.ny + j] = Hzx_yn.data[(i - pis) * Hzx_yn.ny + j] + Hzy_yn.data[i * Hzy_yn.ny + j];
        }
    }
    for (i = pis; i < pie; i++) {
        for (j = 0; j < PMLCellyp; j++) {
            hz.data[i * hz.ny + j + pje] = Hzx_yp.data[(i - pis) * Hzx_yp.ny + j] + Hzy_yp.data[i * Hzy_yp.ny + j];
        }
    }
}

/*******************************************************************************************************
UpdMEltFldForPML_TMz.h
update electric fields at PML regions
TEz
 ***********************************************************************************************************/
void UpdEltFldForPML_TMz(MyStruct ex, MyStruct ey, MyStruct hz) {
    int i, j;
    if (IsPMLxn) {
        for (i = 1; i <= pis; i++)
            for (j = 0; j < ey.ny; j++) {
                ey.data[i * ey.ny + j] = Ceye_xn.data[(i - 1) * Ceye_xn.ny + j] * ey.data[i * ey.ny + j]
                        + Ceyhz_xn.data[(i - 1) * Ceyhz_xn.ny + j]*(hz.data[i * hz.ny + j] - hz.data[(i - 1) * hz.ny + j]);
            }
    }
    if (IsPMLxp) {
        for (i = pie; i < nx; i++)
            for (j = 0; j < ny; j++) {
                ey.data[i * ey.ny + j] = Ceye_xp.data[(i - pie) * Ceye_xp.ny + j] * ey.data[i * ey.ny + j]
                        + Ceyhz_xp.data[(i - pie) * Ceyhz_xp.ny + j]*(hz.data[i * hz.ny + j] - hz.data[(i - 1) * hz.ny + j]);

            }
    }
    if (IsPMLyn) {
        for (i = 0; i <= nxm1; i++)
            for (j = 1; j <= pjs; j++) {
                ex.data[i * ex.ny + j] = Cexe_yn.data[i * Cexe_yn.ny + j - 1] * ex.data[i * ex.ny + j]
                        + Cexhz_yn.data[i * Cexhz_yn.ny + j - 1]*(hz.data[i * hz.ny + j] - hz.data[i * hz.ny + j - 1]);

            }
    }
    if (IsPMLyp) {
        for (i = 0; i <= nxm1; i++)
            for (j = pje; j < ny; j++) {
                ex.data[i * ex.ny + j] = Cexe_yp.data[i * Cexe_yp.ny + j - pje] * ex.data[i * ex.ny + j]
                        + Cexhz_yp.data[i * Cexhz_yp.ny + j - pje]*(hz.data[i * hz.ny + j] - hz.data[i * hz.ny + j - 1]);
            }
    }

}

/*********************************************************************
UpdMEltFldForPML_TEz.h
update electric fields at PML regions
TMz
 *********************************************************************/
void UpdEltFldForPML_TEz(MyStruct ez, MyStruct hx, MyStruct hy) {
    int i, j;
    if (IsPMLxn) {
        for (i = 0; i < Ezx_xn.nx; ++i)
            for (j = 0; j < Ezx_xn.ny; j++) {
                Ezx_xn.data[i * Ezx_xn.ny + j] = Cezxe_xn.data[i * Cezxe_xn.ny + j] * Ezx_xn.data[i * Ezx_xn.ny + j]
                        + Cezxhy_xn.data[i * Cezxhy_xn.ny + j]*(hy.data[(i + 1) * hy.ny + j + 1] - hy.data[i * hy.ny + j + 1]);
            }
        for (i = 0; i < Ezy_xn.nx; ++i)
            for (j = 0; j < Ezy_xn.ny; j++)
                Ezy_xn.data[i * Ezy_xn.ny + j] = Cezye_xn.data[i * Cezye_xn.ny + j] * Ezy_xn.data[i * Ezy_xn.ny + j]
                    + Cezyhx_xn.data[i * Cezyhx_xn.ny + j]*(hx.data[(i + 1) * hx.ny + j + pjs + 1] - hx.data[(i + 1) * hx.ny + j + pjs]);
    }
    /*===============================================================================================*/
    if (IsPMLxp) {
        for (i = 0; i < Ezx_xp.nx; ++i)
            for (j = 0; j < Ezx_xp.ny; j++) {
                Ezx_xp.data[i * Ezx_xp.ny + j] = Cezxe_xp.data[i * Cezxe_xp.ny + j] * Ezx_xp.data[i * Ezx_xp.ny + j]
                        + Cezxhy_xp.data[i * Cezxhy_xp.ny + j]*(hy.data[(i + pie) * hy.ny + j + 1] - hy.data[(i + pie - 1) * hy.ny + j + 1]);
            }
        for (i = 0; i < Ezy_xp.nx; ++i)
            for (j = 0; j < Ezy_xp.ny; j++) {
                Ezy_xp.data[i * Ezy_xp.ny + j] = Cezye_xp.data[i * Cezye_xp.ny + j] * Ezy_xp.data[i * Ezy_xp.ny + j]
                        + Cezyhx_xp.data[i * Cezyhx_xp.ny + j]*(hx.data[(pie + i) * hx.ny + j + pjs + 1] - hx.data[(pie + i) * hx.ny + j + pjs]);
            }
    }
    /*===============================================================================================*/
    if (IsPMLyn) {
        for (i = 0; i < Ezx_yn.nx; ++i)
            for (j = 0; j < Ezx_yn.ny; j++) {
                Ezx_yn.data[i * Ezx_yn.ny + j] = Cezxe_yn.data[i * Cezxe_yn.ny + j] * Ezx_yn.data[i * Ezx_yn.ny + j]
                        + Cezxhy_yn.data[i * Cezxhy_yn.ny + j]*(hy.data[(i + 1 + pis) * hy.ny + j + 1] - hy.data[(i + pis) * hy.ny + j + 1]);
            }
        for (i = 0; i < Ezy_yn.nx; ++i)
            for (j = 0; j < Ezy_yn.ny; j++) {
                Ezy_yn.data[i * Ezy_yn.ny + j] = Cezye_yn.data[i * Cezye_yn.ny + j] * Ezy_yn.data[i * Ezy_yn.ny + j]
                        + Cezyhx_yn.data[i * Cezyhx_yn.ny + j]*(hx.data[(i + 1) * hx.ny + j + 1] - hx.data[(i + 1) * hx.ny + j]);
                /*******************************************************************************************/
            }

    }
    /*===============================================================================================*/
    if (IsPMLyp) {
        for (i = 0; i < Ezx_yp.nx; ++i)
            for (j = 0; j < Ezx_yp.ny; j++) {
                Ezx_yp.data[i * Ezx_yp.ny + j] = Cezxe_yp.data[i * Cezxe_yp.ny + j] * Ezx_yp.data[i * Ezx_yp.ny + j]
                        + Cezxhy_yp.data[i * Cezxhy_yp.ny + j]*(hy.data[(i + pis + 1) * hy.ny + j + pje] - hy.data[(i + pis) * hy.ny + j + pje]);
                /******************************************************************************************/
            }
        for (i = 0; i < Ezy_yp.nx; ++i)
            for (j = 0; j < Ezy_yp.ny; j++) {
                Ezy_yp.data[i * Ezy_yp.ny + j] = Cezye_yp.data[i * Cezye_yp.ny + j] * Ezy_yp.data[i * Ezy_yp.ny + j]
                        + Cezyhx_yp.data[i * Cezyhx_yp.ny + j]*(hx.data[(1 + i) * hx.ny + j + pje] - hx.data[(1 + i) * hx.ny + j + pje - 1]);
                /*******************************************************************************************/
            }
    }
    /*===============================================================================================*/
    for (i = 0; i < pis; i++) {
        for (j = 0; j < pjs; j++) {
            ez.data[(i + 1) * ez.ny + j + 1] = Ezx_xn.data[i * Ezx_xn.ny + j] + Ezy_yn.data[i * Ezy_yn.ny + j];
            /*********************************************************************************************/
        }
    }
    for (i = 0; i < pis; i++) {
        for (j = 0; j < PMLCellyp; j++) {
            ez.data[(i + 1) * ez.ny + j + pje] = Ezx_xn.data[i * Ezx_xn.ny + j + pje - 1] + Ezy_yp.data[i * Ezy_yp.ny + j];
            /*********************************************************************************************/
        }
    }
    for (i = pie; i < nx; i++) {
        for (j = pje; j < ny; j++) {
            ez.data[i * ez.ny + j] = Ezx_xp.data[(i - pie) * Ezx_xp.ny + j - 1] + Ezy_yp.data[(i - 1) * Ezy_yp.ny + j - pje];
            /*********************************************************************************************/
        }
    }
    for (i = pie; i < nx; i++) {
        for (j = 0; j < pjs; j++) {
            ez.data[i * ez.ny + j + 1] = Ezx_xp.data[(i - pie) * Ezx_xp.ny + j] + Ezy_yn.data[(i - 1) * Ezy_yn.ny + j];
            /*********************************************************************************************/
        }
    }
    for (i = pis; i < pie - 1; i++) {
        for (j = 0; j < pjs; j++) {
            ez.data[(i + 1) * ez.ny + j + 1] = Ezx_yn.data[(i - pis) * Ezx_yn.ny + j] + Ezy_yn.data[i * Ezy_yn.ny + j];
            /*********************************************************************************************/
        }
    }
    for (i = pis; i < pie - 1; i++) {
        for (j = 0; j < PMLCellyp; j++) {
            ez.data[(i + 1) * ez.ny + j + pje] = Ezx_yp.data[(i - pis) * Ezx_yp.ny + j] + Ezy_yp.data[i * Ezy_yp.ny + j];
            /*******************************************************************************************************************/
        }
    }
    /*=======================================================================================================*/
    for (i = 0; i < pis; i++) {
        for (j = pjs; j < pje - 1; j++) {
            ez.data[(i + 1) * ez.ny + j + 1] = Ezx_xn.data[i * Ezx_xn.ny + j] + Ezy_xn.data[i * Ezy_xn.ny + j - pjs];
            /*******************************************************************************************************************/
        }
    }
    for (i = 0; i < PMLCellxp; i++) {
        for (j = pjs; j < pje - 1; j++) {
            ez.data[(i + pie) * ez.ny + j + 1] = Ezx_xp.data[i * Ezx_xp.ny + j] + Ezy_xp.data[i * Ezy_xp.ny + j - pjs];
        }
    }
}

void UpdateMFieldForPML(MyStruct hx, MyStruct hy, MyStruct hz, MyStruct ex, MyStruct ey, MyStruct ez) {
    /* update magnetic fields at PML regions */

    if (IsTEx) {
        UpdMagFldForPML_TEz(hx, hy, ez);
    }
    if (IsTMx) {
        UpdMagFldForPML_TMz(hz, ex, ey);
    }

}

void UpdateEFieldForPML(MyStruct hx, MyStruct hy, MyStruct hz, MyStruct ex, MyStruct ey, MyStruct ez) {
    /* update electric fields at PML regions */

    if (IsTEx) {
        UpdEltFldForPML_TEz(ez, hx, hy);
    }
    if (IsTMx) {
        UpdEltFldForPML_TMz(ex, ey, hz);
    }

}

void InitPMLBoundTEx() {
    int i, j, ind;

    MyDataF *rho_e;
    MyDataF *rho_m;

    MyStruct sigma_pey_yn;
    MyStruct sigma_pmy_yn;

    MyStruct sigma_pex_xp;
    MyStruct sigma_pmx_xp;

    MyStruct sigma_pex_xn;
    MyStruct sigma_pmx_xn;

    MyStruct sigma_pey_yp;
    MyStruct sigma_pmy_yp;

    MyDataF sigma_max;

    printf("Initializing pml boundary conditions 2d TEx...\n");
    //printf("nxp1=%d\nnym1=%d\nny=%d\n",nxp1,nym1,ny);
    initMyStruct(&Ezx_xn, PMLCellxn, nym1, 0.0);
    initMyStruct(&Ezx_xp, PMLCellxp, nym1, 0.0);
    initMyStruct(&Ezx_yn, nxm1 - PMLCellxn - PMLCellxp, PMLCellyn, 0.0);
    initMyStruct(&Ezx_yp, nxm1 - PMLCellxn - PMLCellxp, PMLCellyp, 0.0);
    initMyStruct(&Ezy_xn, PMLCellxn, nym1 - PMLCellyn - PMLCellyp, 0.0);
    initMyStruct(&Ezy_xp, PMLCellxp, nym1 - PMLCellyn - PMLCellyp, 0.0);
    initMyStruct(&Ezy_yn, nxm1, PMLCellyn, 0.0);
    initMyStruct(&Ezy_yp, nxm1, PMLCellyp, 0.0);

    if (IsPMLxn) {

        rho_e = (MyDataF*) malloc(PMLCellxn * sizeof (MyDataF));
        rho_m = (MyDataF*) malloc(PMLCellxn * sizeof (MyDataF));

        initMyStruct(&sigma_pex_xn, PMLCellxn, nym1, 0.0);
        initMyStruct(&sigma_pmx_xn, PMLCellxn, nym1, 0.0);

        sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dx * PMLCellxn);

        for (i = PMLCellxn; i > 0; --i) {
            rho_e[PMLCellxn - i] = (i - 0.75) / PMLCellxn;
            rho_m[PMLCellxn - i] = (i - 0.25) / PMLCellxn;
        }

        for (i = 0; i < PMLCellxn; i++)
            for (j = 0; j < nym1; ++j) {
                sigma_pex_xn.data[i * nym1 + j] = sigma_max * pow(rho_e[i], pml_order);
                sigma_pmx_xn.data[i * nym1 + j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[i], pml_order);
            }


        initMyStruct(&Chyh_xn, PMLCellxn, nym1, 0.0);
        initMyStruct(&Chyez_xn, PMLCellxn, nym1, 0.0);

        initMyStruct(&Cezxe_xn, Ezx_xn.nx, Ezx_xn.ny, 0.0);
        initMyStruct(&Cezxhy_xn, Ezx_xn.nx, Ezx_xn.ny, 0.0);

        initMyStruct(&Cezye_xn, Ezy_xn.nx, Ezy_xn.ny, 1.0);
        initMyStruct(&Cezyhx_xn, Ezy_xn.nx, Ezy_xn.ny, -dt / (dy * eps_0));

        ind = 0;
        for (i = 0; i < PMLCellxn; ++i)
            for (j = 0; j < nym1; j++) {
                /* Coefficients updating Hx */
                Chyh_xn.data[ind] = (2 * mu_0 - dt * sigma_pmx_xn.data[ind]) / (2 * mu_0 + dt * sigma_pmx_xn.data[ind]);
                Chyez_xn.data[ind] = (2 * dt / dx) / (2 * mu_0 + dt * sigma_pmx_xn.data[ind]);

                /* Coefficients updating Ezx */
                Cezxe_xn.data[ind] = (2 * eps_0 - dt * sigma_pex_xn.data[ind]) / (2 * eps_0 + dt * sigma_pex_xn.data[ind]);
                Cezxhy_xn.data[ind] = (2 * dt / dx) / (2 * eps_0 + dt * sigma_pex_xn.data[ind]);

                ind++;
            }

        free(sigma_pmx_xn.data);
        free(sigma_pex_xn.data);

        free(rho_e);
        free(rho_m);
    }

    if (IsPMLyn) {

        rho_e = (MyDataF*) malloc(PMLCellyn * sizeof (MyDataF));
        rho_m = (MyDataF*) malloc(PMLCellyn * sizeof (MyDataF));

        initMyStruct(&sigma_pey_yn, nxm1, PMLCellyn, 0.0);
        initMyStruct(&sigma_pmy_yn, nxm1, PMLCellyn, 0.0);

        sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dy * PMLCellyn);

        for (i = PMLCellyn; i > 0; --i) {
            rho_e[PMLCellyn - i] = (i - 0.75) / PMLCellyn;
            rho_m[PMLCellyn - i] = (i - 0.25) / PMLCellyn;
        }

        for (i = 0; i < PMLCellyn; i++)
            for (j = 0; j < nxm1; ++j) {
                sigma_pey_yn.data[j * PMLCellyn + i] = sigma_max * pow(rho_e[i], pml_order);
                sigma_pmy_yn.data[j * PMLCellyn + i] = (mu_0 / eps_0) * sigma_max * pow(rho_m[i], pml_order);
            }

        initMyStruct(&Chxh_yn, nxm1, PMLCellyn, 0.0);
        initMyStruct(&Chxez_yn, nxm1, PMLCellyn, 0.0);


        initMyStruct(&Cezye_yn, Ezy_yn.nx, Ezy_yn.ny, 0.0);
        initMyStruct(&Cezyhx_yn, Ezy_yn.nx, Ezy_yn.ny, 0.0);

        ind = 0;
        for (i = 0; i < nxm1; ++i)
            for (j = 0; j < PMLCellyn; j++) {
                /* Coefficients updating Hx */
                Chxh_yn.data[ind] = (2 * mu_0 - dt * sigma_pmy_yn.data[ind]) / (2 * mu_0 + dt * sigma_pmy_yn.data[ind]);
                Chxez_yn.data[ind] = -(2 * dt / dy) / (2 * mu_0 + dt * sigma_pmy_yn.data[ind]);

                /* Coefficients updating Ezy */
                Cezye_yn.data[ind] = (2 * eps_0 - dt * sigma_pey_yn.data[ind]) / (2 * eps_0 + dt * sigma_pey_yn.data[ind]);
                Cezyhx_yn.data[ind] = -(2 * dt / dy) / (2 * eps_0 + dt * sigma_pey_yn.data[ind]);
                ind++;
            }


        initMyStruct(&Cezxe_yn, Ezx_yn.nx, Ezx_yn.ny, 1.0);
        initMyStruct(&Cezxhy_yn, Ezx_yn.nx, Ezx_yn.ny, dt / (dx * eps_0));

        free(sigma_pmy_yn.data);
        free(sigma_pey_yn.data);

        free(rho_e);
        free(rho_m);
    }

    if (IsPMLxp) {

        rho_e = (MyDataF*) malloc(PMLCellxp * sizeof (MyDataF));
        rho_m = (MyDataF*) malloc(PMLCellxp * sizeof (MyDataF));

        initMyStruct(&sigma_pex_xp, PMLCellxp, nym1, 0.0);
        initMyStruct(&sigma_pmx_xp, PMLCellxp, nym1, 0.0);

        sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dx * PMLCellxp);

        for (i = 0; i < PMLCellyp; ++i) {
            rho_e[i] = (i + 1 - 0.75) / PMLCellxp;
            rho_m[i] = (i + 1 - 0.25) / PMLCellxp;
        }

        for (i = 0; i < nym1; i++)
            for (j = 0; j < PMLCellxp; ++j) {
                sigma_pex_xp.data[j * nym1 + i] = sigma_max * pow(rho_e[j], pml_order);
                sigma_pmx_xp.data[j * nym1 + i] = (mu_0 / eps_0) * sigma_max * pow(rho_m[j], pml_order);
            }

        initMyStruct(&Chyh_xp, PMLCellxp, nym1, 0.0);
        initMyStruct(&Chyez_xp, PMLCellxp, nym1, 0.0);

        initMyStruct(&Cezxe_xp, Ezx_xp.nx, Ezx_xp.ny, 0.0);
        initMyStruct(&Cezxhy_xp, Ezx_xp.nx, Ezx_xp.ny, 0.0);

        ind = 0;
        for (i = 0; i < PMLCellxp; ++i)
            for (j = 0; j < nym1; j++) {
                /* Coefficients updating Hy */
                Chyh_xp.data[ind] = (2 * mu_0 - dt * sigma_pmx_xp.data[ind]) / (2 * mu_0 + dt * sigma_pmx_xp.data[ind]);
                Chyez_xp.data[ind] = (2 * dt / dx) / (2 * mu_0 + dt * sigma_pmx_xp.data[ind]);

                /* Coefficients updating Ezx */
                Cezxe_xp.data[ind] = (2 * eps_0 - dt * sigma_pex_xp.data[ind]) / (2 * eps_0 + dt * sigma_pex_xp.data[ind]);
                Cezxhy_xp.data[ind] = (2 * dt / dx) / (2 * eps_0 + dt * sigma_pex_xp.data[ind]);

                ind++;
            }


        initMyStruct(&Cezye_xp, Ezy_xp.nx, Ezy_xp.ny, 1.0);
        initMyStruct(&Cezyhx_xp, Ezy_xp.nx, Ezy_xp.ny, -dt / (dy * eps_0));

        free(sigma_pmx_xp.data);
        free(sigma_pex_xp.data);

        free(rho_e);
        free(rho_m);
    }
    if (IsPMLyp) {

        rho_e = (MyDataF*) malloc(PMLCellyp * sizeof (MyDataF));
        rho_m = (MyDataF*) malloc(PMLCellyp * sizeof (MyDataF));

        initMyStruct(&sigma_pey_yp, nxm1, PMLCellyp, 0.0);
        initMyStruct(&sigma_pmy_yp, nxm1, PMLCellyp, 0.0);

        sigma_max = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / (2 * dy * PMLCellyp);

        for (i = 0; i < PMLCellyp; ++i) {
            rho_e[i] = (i + 1 - 0.75) / PMLCellyp;
            rho_m[i] = (i + 1 - 0.25) / PMLCellyp;
        }

        for (i = 0; i < nxm1; i++)
            for (j = 0; j < PMLCellyp; ++j) {
                sigma_pey_yp.data[i * PMLCellyp + j] = sigma_max * pow(rho_e[j], pml_order);
                sigma_pmy_yp.data[i * PMLCellyp + j] = (mu_0 / eps_0) * sigma_max * pow(rho_m[j], pml_order);
            }

        initMyStruct(&Chxh_yp, nxm1, PMLCellyp, 0.0);
        initMyStruct(&Chxez_yp, nxm1, PMLCellyp, 0.0);

        initMyStruct(&Cezye_yp, Ezy_yp.nx, Ezy_yp.ny, 0.0);
        initMyStruct(&Cezyhx_yp, Ezy_yp.nx, Ezy_yp.ny, 0.0);

        ind = 0;
        for (i = 0; i < nxm1; ++i)
            for (j = 0; j < PMLCellyp; j++) {
                /* Coefficients updating Hy */
                Chxh_yp.data[ind] = (2 * mu_0 - dt * sigma_pmy_yp.data[ind]) / (2 * mu_0 + dt * sigma_pmy_yp.data[ind]);
                Chxez_yp.data[ind] = -(2 * dt / dx) / (2 * mu_0 + dt * sigma_pmy_yp.data[ind]);

                /* Coefficients updating Ezy */
                Cezye_yp.data[ind] = (2 * eps_0 - dt * sigma_pey_yp.data[ind]) / (2 * eps_0 + dt * sigma_pey_yp.data[ind]);
                Cezyhx_yp.data[ind] = -(2 * dt / dy) / (2 * eps_0 + dt * sigma_pey_yp.data[ind]);
                ind++;
            }

        initMyStruct(&Cezxe_yp, Ezx_yp.nx, Ezx_yp.ny, 1.0);
        initMyStruct(&Cezxhy_yp, Ezx_yp.nx, Ezx_yp.ny, dt / (dx * eps_0));

        free(sigma_pmy_yp.data);
        free(sigma_pey_yp.data);

        free(rho_e);
        free(rho_m);
    }


    printf("End of initializing pml boundary conditions 2d TEx...\n");
}

void InitPMLBoundTMx() {

    int i, j, ind;
    MyDataF pmx, pex;
    MyDataF rho_e;
    MyDataF rho_m;
    MyDataF mu_d_eps = mu_0 / eps_0;
    MyDataF sigma_max;
    MyDataF tsigmax = -(pml_order + 1) * eps_0 * c * log(pml_R_0) / 2;
    MyDataF dyepsdt = dy * eps_0 / dt;
    MyDataF dxepsdt = dx * eps_0 / dt;
    MyDataF dxmudt = dx * mu_0 / dt;
    MyDataF dymudt = dy * mu_0 / dt;
    MyDataF hdx = dx / 2;
    MyDataF hdy = dy / 2;
    printf("Initializing pml boundary conditions 2d TMx...\n");

    initMyStruct(&Hzx_xn, PMLCellxn, ny, 0.0);
    initMyStruct(&Hzx_xp, PMLCellxp, ny, 0.0);

    initMyStruct(&Hzx_yn, nx - PMLCellxn - PMLCellxp, PMLCellyn, 0.0);
    initMyStruct(&Hzx_yp, nx - PMLCellxp - PMLCellxn, PMLCellyp, 0.0);

    initMyStruct(&Hzy_xn, PMLCellxn, ny - PMLCellyn - PMLCellyp, 0.0);
    initMyStruct(&Hzy_xp, PMLCellxp, ny - PMLCellyp - PMLCellyn, 0.0);

    initMyStruct(&Hzy_yn, nx, PMLCellyn, 0.0);
    initMyStruct(&Hzy_yp, nx, PMLCellyp, 0.0);

    /*===============================================================================================*/
    if (IsPMLxn) {

        initMyStruct(&Ceye_xn, PMLCellxn, ny, 0.0);
        initMyStruct(&Ceyhz_xn, PMLCellxn, ny, 0.0);

        initMyStruct(&Chzxh_xn, Hzx_xn.nx, Hzx_xn.ny, 0.0);
        initMyStruct(&Chzxey_xn, Hzx_xn.nx, Hzx_xn.ny, 0.0);

        initMyStruct(&Chzyh_xn, Hzy_xn.nx, Hzy_xn.ny, 1.0);
        initMyStruct(&Chzyex_xn, Hzy_xn.nx, Hzy_xn.ny, dt / (dy * mu_0));

        sigma_max = tsigmax / (dx * PMLCellxn);

        for (ind = i = 0; i < PMLCellxn; ++i) {
            rho_e = 1 - (i - 0.75) / PMLCellxn;
            rho_m = 1 - (i - 0.25) / PMLCellxn;
            pex = sigma_max * pow(rho_e, pml_order);
            pmx = mu_d_eps * sigma_max * pow(rho_m, pml_order);
            for (j = 0; j < ny; j++) {
                /* Coefficients updating Ey */
                Ceye_xn.data[ind] = (2 * eps_0 - dt * pex) / (2 * eps_0 + dt * pex);
                Ceyhz_xn.data[ind] = -1 / (dxepsdt + hdx * pex);

                /* Coefficients updating Hzx */
                Chzxh_xn.data[ind] = (2 * mu_0 - dt * pmx) / (2 * mu_0 + dt * pmx);
                Chzxey_xn.data[ind++] = -1 / (dxmudt + hdx * pmx);
            }
        }
    }

    /*===============================================================================================*/
    if (IsPMLxp) {

        sigma_max = tsigmax / (dx * PMLCellxp);

        initMyStruct(&Ceye_xp, PMLCellxp, ny, 0.0);
        initMyStruct(&Ceyhz_xp, PMLCellxp, ny, 0.0);

        initMyStruct(&Chzxh_xp, Hzx_xp.nx, Hzx_xp.ny, 0.0);
        initMyStruct(&Chzxey_xp, Hzx_xp.nx, Hzx_xp.ny, 0.0);

        initMyStruct(&Chzyh_xp, Hzy_xp.nx, Hzy_xp.ny, 1.0);
        initMyStruct(&Chzyex_xp, Hzy_xp.nx, Hzy_xp.ny, dt / (dy * mu_0));

        for (ind = i = 0; i < PMLCellxp; ++i) {
            rho_e = (i + 0.25) / PMLCellxp;
            rho_m = (i + 0.75) / PMLCellxp;
            pex = sigma_max * pow(rho_e, pml_order);
            pmx = mu_d_eps * sigma_max * pow(rho_m, pml_order);

            for (j = 0; j < ny; j++) {
                /* Coefficients updating Ey */
                Ceye_xp.data[ind] = (2 * eps_0 - dt * pex) / (2 * eps_0 + dt * pex);
                Ceyhz_xp.data[ind] = -1 / (dxepsdt + hdx * pex);

                /* Coefficients updating Hzx */
                Chzxh_xp.data[ind] = (2 * mu_0 - dt * pmx) / (2 * mu_0 + dt * pmx);
                Chzxey_xp.data[ind] = -1 / (dxmudt + hdx * pmx);
                ind++;
            }
        }
    }
    /*================================================================================================*/
    if (IsPMLyn) {

        sigma_max = tsigmax / (dy * PMLCellyn);


        initMyStruct(&Cexe_yn, nx, PMLCellyn, 0.0);
        initMyStruct(&Cexhz_yn, nx, PMLCellyn, 0.0);

        initMyStruct(&Chzxh_yn, Hzx_yn.nx, Hzx_yn.ny, 1.0);
        initMyStruct(&Chzxey_yn, Hzx_yn.nx, Hzx_yn.ny, -dt / (dx * mu_0));

        initMyStruct(&Chzyh_yn, Hzy_yn.nx, Hzy_yn.ny, 0.0);
        initMyStruct(&Chzyex_yn, Hzy_yn.nx, Hzy_yn.ny, 0.0);


        for (j = 0; j < PMLCellyn; j++) {
            rho_e = 1 - (j - 0.75) / PMLCellyn;
            rho_m = 1 - (j - 0.25) / PMLCellyn;
            pex = sigma_max * pow(rho_e, pml_order);
            pmx = mu_d_eps * sigma_max * pow(rho_m, pml_order);
            ind = j;
            for (i = 0; i < nx; ++i) {
                /* Coefficients updating Ey */
                Cexe_yn.data[ind] = (2 * eps_0 - dt * pex) / (2 * eps_0 + dt * pex);
                Cexhz_yn.data[ind] = 1 / (dyepsdt + hdy * pex);

                /* Coefficients updating Hzy */
                Chzyh_yn.data[ind] = (2 * mu_0 - dt * pmx) / (2 * mu_0 + dt * pmx);
                Chzyex_yn.data[ind] = 1 / (dymudt + hdy * pmx);
                ind += PMLCellyn;
            }
        }
    }
    /*===============================================================================================*/
    if (IsPMLyp) {

        sigma_max = tsigmax / (dy * PMLCellyp);

        initMyStruct(&Cexe_yp, nx, PMLCellyp, 0.0);
        initMyStruct(&Cexhz_yp, nx, PMLCellyp, 0.0);

        initMyStruct(&Chzxh_yp, Hzx_yp.nx, Hzx_yp.ny, 1.0);
        initMyStruct(&Chzxey_yp, Hzx_yp.nx, Hzx_yp.ny, -dt / (dx * mu_0));

        initMyStruct(&Chzyh_yp, Hzy_yp.nx, Hzy_yp.ny, 0.0);
        initMyStruct(&Chzyex_yp, Hzy_yp.nx, Hzy_yp.ny, 0.0);

        for (j = 0; j < PMLCellyp; j++) {
            rho_e = (j + 0.25) / PMLCellyp;
            rho_m = (j + 0.75) / PMLCellyp;
            pex = sigma_max * pow(rho_e, pml_order);
            pmx = mu_d_eps * sigma_max * pow(rho_m, pml_order);
            ind = j;
            for (i = 0; i < nx; ++i) {
                /* Coefficients updating Ey */
                Cexe_yp.data[ind] = (2 * eps_0 - dt * pex) / (2 * eps_0 + dt * pex);
                Cexhz_yp.data[ind] = 1 / (dyepsdt + hdy * pex);

                /* Coefficients updating Hzy */
                Chzyh_yp.data[ind] = (2 * mu_0 - dt * pmx) / (2 * mu_0 + dt * pmx);
                Chzyex_yp.data[ind] = 1 / (dymudt + hdy * pmx);
                ind += PMLCellyp;
            }
        }
    }

    printf("End of initializing pml boundary conditions 2d TMx...\n");
}

void InitPMLBound() {

    /* determine the boundaries of the non-pml region */
    if (IsAnySidePML) {
        if (IsTMx)
            InitPMLBoundTMx();
        if (IsTEx)
            InitPMLBoundTEx();
    }
}

void FreePMLSpace() {
    if (IsTMx) {

        printf("Freeing TMx PML space...\n");

        free(Hzx_xn.data);
        free(Hzx_xp.data);

        free(Hzx_yn.data);
        free(Hzx_yp.data);

        free(Hzy_xn.data);
        free(Hzy_xp.data);

        free(Hzy_yn.data);
        free(Hzy_yp.data);

        /*=====================================================================================================================*/
        if (IsPMLxn) {
            free(Ceye_xn.data);
            free(Ceyhz_xn.data);

            free(Chzxh_xn.data);
            free(Chzxey_xn.data);

            free(Chzyh_xn.data);
            free(Chzyex_xn.data);
        }

        /*====================================================================================================================*/
        if (IsPMLxp) {
            free(Ceye_xp.data);
            free(Ceyhz_xp.data);

            free(Chzxh_xp.data);
            free(Chzxey_xp.data);

            free(Chzyh_xp.data);
            free(Chzyex_xp.data);
        }
        /*=====================================================================================================================*/
        if (IsPMLyn) {

            free(Cexe_yn.data);
            free(Cexhz_yn.data);

            free(Chzxh_yn.data);
            free(Chzxey_yn.data);

            free(Chzyh_yn.data);
            free(Chzyex_yn.data);


        }
        /*=====================================================================================================================*/
        if (IsPMLyp) {


            free(Cexe_yp.data);
            free(Cexhz_yp.data);

            free(Chzxh_yp.data);
            free(Chzxey_yp.data);

            free(Chzyh_yp.data);
            free(Chzyex_yp.data);

        }

        printf("End of freeing TMx  PML space...\n");
    }

    if (IsTEx) {


        printf("Freeing TEx PML space...\n");

        free(Ezx_xn.data);
        free(Ezx_xp.data);
        free(Ezx_yn.data);
        free(Ezx_yp.data);
        free(Ezy_xn.data);
        free(Ezy_xp.data);
        free(Ezy_yn.data);
        free(Ezy_yp.data);

        if (IsPMLxn) {
            free(Chyh_xn.data);
            free(Chyez_xn.data);

            free(Cezxe_xn.data);
            free(Cezxhy_xn.data);

            free(Cezye_xn.data);
            free(Cezyhx_xn.data);
        }

        if (IsPMLyn) {


            free(Chxh_yn.data);
            free(Chxez_yn.data);


            free(Cezye_yn.data);
            free(Cezyhx_yn.data);


            free(Cezxe_yn.data);
            free(Cezxhy_yn.data);

        }

        if (IsPMLxp) {
            free(Chyh_xp.data);
            free(Chyez_xp.data);

            free(Cezxe_xp.data);
            free(Cezxhy_xp.data);


            free(Cezye_xp.data);
            free(Cezyhx_xp.data);
        }
        if (IsPMLyp) {
            free(Chxh_yp.data);
            free(Chxez_yp.data);

            free(Cezye_yp.data);
            free(Cezyhx_yp.data);

            free(Cezxe_yp.data);
            free(Cezxhy_yp.data);

        }
        printf("End of freeing TEx PML Space...\n");
    }
}