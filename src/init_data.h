

#ifndef INIT_DATA_INCLUDE
#define INIT_DATA_INCLUDE

//#include"init_com_data.h"
int cnepx, cnepy;

void Init_ne() {
    MyDataF temp = 2500e-12;
    int i, j, ind;
    const int bndsz = NUMBER_OF_CELLS_IN_PML_BOUND + SCATTER_FIELD_DOMAIN_BND_SIZE;
    MyDataF dx2, dy2;
    dx2 = dx / m;
    dy2 = dy / m;

    for (i = 0; i < ne.nx; i++)
        for (j = 0; j < ne.ny; j++) {
            ind = i * ne.ny + j;
            ne_pre.data[ind] = ne.data[ind] = //1e6;
                    NE0 * exp(-(((i - bndsz * m) * dx2 - 2.25 * lamda)*((i - bndsz * m) * dx2 - 2.25 * lamda)+(j - 0.5 * ne.ny)*(j - 0.5 * ne.ny) * dy2 * dy2) / temp);

        }
}

void InitBeta(MyStruct* stru, MyDataF gamma) {
    int i, j, ind;
    MyDataF temp = 0.25 * e * e * dt * dt / me / eps_0 / gamma;
    for (i = 0; i < stru->nx; i++)
        for (j = 0; j < stru->ny; j++) {
            ind = i * stru->ny + j;
            stru->data[ind] = temp * ne.data[ind];
        }
}

void InitCee(MyStruct beta) {
    int i, j, index;
    if (IsTMx) {
        for (i = 0; i < Ceex.nx; ++i)
            for (j = 0; j < Ceex.ny; j++) {
                index = (m * i + m / 2) * beta.ny + m*j;
                Ceex.data[i * Ceex.ny + j] = -beta.data[index] / (1 + beta.data[index]) + 1 / (1 + beta.data[index]);
            }
        for (i = 0; i < Ceey.ny; ++i)
            for (j = 0; j < Ceey.ny; j++) {
                index = m * i * beta.ny + m * j + m / 2;
                Ceey.data[i * Ceey.ny + j] = -beta.data[index] / (1 + beta.data[index]) + 1 / (1 + beta.data[index]);
            }
    }
    if (IsTEx) {
        for (i = 0; i < Ceez.nx; ++i)
            for (j = 0; j < Ceez.ny; j++) {
                index = m * i * beta.ny + m*j;
                Ceez.data[i * Ceez.ny + j] = -beta.data[index] / (1 + beta.data[index]) + 1 / (1 + beta.data[index]);
            }
    }
    //PrintData(Ceez);
}

void InitCeh(MyStruct beta) {
    int i, j, index = 0;
    //MyDataF temp=dt/eps_0;
    if (IsTMx) {
        for (i = 0; i < Cehx.nx; ++i)
            for (j = 0; j < Cehx.ny; j++) {
                index = (m * i + m / 2) * beta.ny + m*j;
                Cehx.data[i * Cehx.ny + j] = dt / (eps_0 + eps_0 * beta.data[index]) / dy;
            }
        for (i = 0; i < Cehy.ny; ++i)
            for (j = 0; j < Cehy.ny; j++) {
                index = m * i * beta.ny + m * j + m / 2;
                Cehy.data[i * Cehy.ny + j] = -dt / (eps_0 + eps_0 * beta.data[index]) / dx;
            }
    }
    if (IsTEx) {
        for (i = 0; i < Cehz.nx; ++i)
            for (j = 0; j < Cehz.ny; j++) {
                index = m * i * beta.ny + m*j;
                Cehz.data[i * Cehz.ny + j] = dt / (eps_0 + eps_0 * beta.data[index]) / dx;
            }
    }
    //PrintData(Cehz);
}

void InitCev(MyStruct beta, MyDataF alpha_ev) {
    int i, j, index = 0;
    MyDataF temp;

    temp = 0.5 * e * dt * (1 + alpha_ev) / eps_0;
    if (IsTMx) {
        for (i = 0; i < Cevx.nx; ++i)
            for (j = 0; j < Cevx.ny; j++) {
                index = (m * i + m / 2) * beta.ny + 2 * j;
                Cevx.data[i * Cevx.ny + j] = temp * ne.data[index] / (1 + beta.data[index]);
            }
        for (i = 0; i < Cevy.ny; ++i)
            for (j = 0; j < Cevy.ny; j++) {
                index = m * i * beta.ny + m * j + m / 2;
                Cevy.data[i * Cevy.ny + j] = temp * ne.data[index] / (1 + beta.data[index]);
            }
    }
    if (IsTEx) {
        for (i = 0; i < Cevz.nx; ++i)
            for (j = 0; j < Cevz.ny; j++) {
                index = m * i * beta.ny + m*j;
                Cevz.data[i * Cevz.ny + j] = temp * ne.data[index] / (1 + beta.data[index]);
            }
        //PrintData(Cevz);
    }
}
/****************************************************************/
//coefficients when without density

void InitCeeNoDen() {
    int i, j, index;
    if (IsTMx) {
        for (i = 0; i < Ceex.nx; ++i)
            for (j = 0; j < Ceex.ny; j++) {
                Ceex.data[i * Ceex.ny + j] = 1;
            }
        for (i = 0; i < Ceey.ny; ++i)
            for (j = 0; j < Ceey.ny; j++) {
                Ceey.data[i * Ceey.ny + j] = 1;
            }
    }
    if (IsTEx) {
        for (i = 0; i < Ceez.nx; ++i)
            for (j = 0; j < Ceez.ny; j++) {
                Ceez.data[i * Ceez.ny + j] = 1;
            }
    }
    //PrintData(Ceez);
}

void InitCehNoDen() {
    int i, j, index = 0;
    //MyDataF temp=dt/eps_0;
    if (IsTMx) {
        for (i = 0; i < Cehx.nx; ++i)
            for (j = 0; j < Cehx.ny; j++) {
                Cehx.data[i * Cehx.ny + j] = dt / (eps_0) / dy;
            }
        for (i = 0; i < Cehy.ny; ++i)
            for (j = 0; j < Cehy.ny; j++) {
                Cehy.data[i * Cehy.ny + j] = -dt / (eps_0) / dx;
            }
    }
    if (IsTEx) {
        for (i = 0; i < Cehz.nx; ++i)
            for (j = 0; j < Cehz.ny; j++) {
                Cehz.data[i * Cehz.ny + j] = dt / (eps_0) / dx;
            }
    }
    //PrintData(Cehz);
}

void InitCevNoDen() {
    int i, j, index = 0;
    MyDataF temp;
    if (IsTMx) {
        for (i = 0; i < Cevx.nx; ++i)
            for (j = 0; j < Cevx.ny; j++) {
                Cevx.data[i * Cevx.ny + j] = 0;
            }
        for (i = 0; i < Cevy.ny; ++i)
            for (j = 0; j < Cevy.ny; j++) {
                Cevy.data[i * Cevy.ny + j] = 0;
            }
    }
    if (IsTEx) {
        for (i = 0; i < Cevz.nx; ++i)
            for (j = 0; j < Cevz.ny; j++) {
                Cevz.data[i * Cevz.ny + j] = 0;
            }
        //PrintData(Cevz);
    }
}

/**********************************************************************************************/

void InitCoeff() {
    MyDataF gamma, a = 0;

    a = 0.5 * dt*vm;
    gamma = 1 + a;
    alpha = (1 - a) / (1 + a);

    Chxez = -dt / mu_0 / dy;
    Chyez = dt / mu_0 / dx;
    Chzex = dt / mu_0 / dy;
    Chzey = -dt / mu_0 / dx;
    Cve = e * dt * 0.5 / (me * gamma);
    if (IfWithDensity) {
        Init_ne();
        InitBeta(&beta, gamma);
        /*PrintData(beta);*/
        InitCee(beta);
        InitCeh(beta);
        InitCev(beta, alpha);
    } else {
        InitCeeNoDen();
        InitCehNoDen();
        InitCevNoDen();
    }
}
//void InitEMField(){
//	if(IsTMx){
//		ResetStructData(Ex);
//		ResetStructData(Ey);
//		ResetStructData(Hz);
//
//		ResetStructData(Ex_s);
//		ResetStructData(Ey_s);
//		ResetStructData(Hz_s);
//		ResetStructData(Hz_s_pre);
//
//		ResetStructData(Ex_i_pre);
//		ResetStructData(Ey_i_pre);
//		ResetStructData(Ex_s_pre);
//		ResetStructData(Ey_s_pre);
//
//		ResetStructData(Ex_i);
//		ResetStructData(Ey_i);
//		ResetStructData(Hz_i);
//		ResetStructData(Hz_i_pre);
//
//		ResetStructData(Vex);
//		ResetStructData(Vey);
//	}
//	if(IsTEx){
//		ResetStructData(Hx);
//		ResetStructData(Hy);
//		ResetStructData(Ez);
//
//		ResetStructData(Hx_s);
//		ResetStructData(Hy_s);
//		ResetStructData(Hx_s_pre);
//		ResetStructData(Hy_s_pre);
//		ResetStructData(Ez_s);
//
//		ResetStructData(Hx_i);
//		ResetStructData(Hy_i);
//		ResetStructData(Ez_i);
//
//		ResetStructData(Hx_i_pre);
//		ResetStructData(Hy_i_pre);
//		ResetStructData(Ez_s_pre);
//		ResetStructData(Ez_i_pre);
//		
//		ResetStructData(Vez);	
//	}
//	ResetStructData(Ez);
//}
//

#endif
