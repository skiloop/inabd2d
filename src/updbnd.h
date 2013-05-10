#ifndef __UPDBND_H__
#define __UPDBND_H__

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
    /*===================================================================================================================*/
    for (i = 0; i < pis; i++) {
        for (j = 0; j < pjs; j++) {
            ez.data[(i + 1) * ez.ny + j + 1] = Ezx_xn.data[i * Ezx_xn.ny + j] + Ezy_yn.data[i * Ezy_yn.ny + j];
            /*******************************************************************************************************************/
        }
    }
    for (i = 0; i < pis; i++) {
        for (j = 0; j < PMLCellyp; j++) {
            ez.data[(i + 1) * ez.ny + j + pje] = Ezx_xn.data[i * Ezx_xn.ny + j + pje - 1] + Ezy_yp.data[i * Ezy_yp.ny + j];
            /*******************************************************************************************************************/
        }
    }
    for (i = pie; i < nx; i++) {
        for (j = pje; j < ny; j++) {
            ez.data[i * ez.ny + j] = Ezx_xp.data[(i - pie) * Ezx_xp.ny + j - 1] + Ezy_yp.data[(i - 1) * Ezy_yp.ny + j - pje];
            /*******************************************************************************************************************/
        }
    }
    for (i = pie; i < nx; i++) {
        for (j = 0; j < pjs; j++) {
            ez.data[i * ez.ny + j + 1] = Ezx_xp.data[(i - pie) * Ezx_xp.ny + j] + Ezy_yn.data[(i - 1) * Ezy_yn.ny + j];
            /*******************************************************************************************************************/
        }
    }
    for (i = pis; i < pie - 1; i++) {
        for (j = 0; j < pjs; j++) {
            ez.data[(i + 1) * ez.ny + j + 1] = Ezx_yn.data[(i - pis) * Ezx_yn.ny + j] + Ezy_yn.data[i * Ezy_yn.ny + j];
            /*******************************************************************************************************************/
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
#endif
