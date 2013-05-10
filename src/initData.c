
#include <math.h>

#include "common.h"
#include "dataType.h"
#include "commonData.h"
#include "fdtd.h"
#include "InonizationFormula.h"

void Init_ne() {
    MyDataF temp = 2500e-12;
    int i, j, ind;
    const int bndsz = NUMBER_OF_CELLS_IN_PML_BOUND + SCATTER_FIELD_DOMAIN_BND_SIZE;
    MyDataF dx2, dy2;
    dx2 = dx / m;
    dy2 = dy / m;

    for (i = 0; i < ne.nx; i++) {
        for (j = 0; j < ne.ny; j++) {
            ind = i * ne.ny + j;
            ne_pre.data[ind] = ne.data[ind] = //1e6;
                    NE0 * exp(-(((i - bndsz * m) * dx2 - 2.25 * lamda)*((i - bndsz * m) * dx2 - 2.25 * lamda)+(j - 0.5 * ne.ny)*(j - 0.5 * ne.ny) * dy2 * dy2) / temp);
        }
    }
}

void InitBeta(MyStruct* stru, MyDataF gamma) {
    int i, j, ind;
    MyDataF temp = 0.25 * e * e * dt * dt / me / eps_0 / gamma, a, vm_, gamma_;
    for (i = 0; i < stru->nx; i++) {
        for (j = 0; j < stru->ny; j++) {
            ind = i * stru->ny + j;
            if (niutype == 4) {
                vm_ = collision(Erms.data[ind], p);
                a = dt * vm_ / 2;
                gamma_ = 1 + a;
                temp = 0.25 * e * e * dt * dt / me / eps_0 / gamma_;
            } else if (niutype == 5) {
                a = dt * Nu_c.data[ind] / 2;
                gamma_ = 1 + a;
                temp = 0.25 * e * e * dt * dt / me / eps_0 / gamma_;
            }
            stru->data[ind] = temp * ne.data[ind];
        }
    }
}

void InitCee(MyStruct beta) {
    int i, j, index;
    if (IsTMx) {
        for (i = 0; i < Ceex.nx; ++i) {
            for (j = 0; j < Ceex.ny; j++) {
                index = (m * i + m / 2) * beta.ny + m*j;
                Ceex.data[i * Ceex.ny + j] = -beta.data[index] / (1 + beta.data[index]) + 1 / (1 + beta.data[index]);
            }
        }
        for (i = 0; i < Ceey.ny; ++i) {
            for (j = 0; j < Ceey.ny; j++) {
                index = m * i * beta.ny + m * j + m / 2;
                Ceey.data[i * Ceey.ny + j] = -beta.data[index] / (1 + beta.data[index]) + 1 / (1 + beta.data[index]);
            }
        }
    }
    if (IsTEx) {
        for (i = 0; i < Ceez.nx; ++i) {
            for (j = 0; j < Ceez.ny; j++) {
                index = m * i * beta.ny + m*j;
                Ceez.data[i * Ceez.ny + j] = -beta.data[index] / (1 + beta.data[index]) + 1 / (1 + beta.data[index]);
            }
        }
    }
    //PrintData(Ceez);
}

void InitCeh(MyStruct beta) {
    int i, j, index = 0;
    //MyDataF temp=dt/eps_0;
    if (IsTMx) {
        for (i = 0; i < Cehx.nx; ++i) {
            for (j = 0; j < Cehx.ny; j++) {
                index = (m * i + m / 2) * beta.ny + m*j;
                Cehx.data[i * Cehx.ny + j] = dt / (eps_0 + eps_0 * beta.data[index]) / dy;
            }
        }
        for (i = 0; i < Cehy.ny; ++i) {
            for (j = 0; j < Cehy.ny; j++) {
                index = m * i * beta.ny + m * j + m / 2;
                Cehy.data[i * Cehy.ny + j] = -dt / (eps_0 + eps_0 * beta.data[index]) / dx;
            }
        }
    }
    if (IsTEx) {
        for (i = 0; i < Cehz.nx; ++i) {
            for (j = 0; j < Cehz.ny; j++) {
                index = m * i * beta.ny + m*j;
                Cehz.data[i * Cehz.ny + j] = dt / (eps_0 + eps_0 * beta.data[index]) / dx;
            }
        }
    }
    //PrintData(Cehz);
}

void InitCev(MyStruct beta, MyDataF alpha_ev) {
    int i, j, index = 0;
    MyDataF temp, a, vm_;

    temp = 0.5 * e * dt * (1 + alpha_ev) / eps_0;
    if (IsTMx) {
        for (i = 0; i < Cevx.nx; ++i) {
            for (j = 0; j < Cevx.ny; j++) {
                index = (m * i + m / 2) * beta.ny + 2 * j;
                if (niutype == 4) {
                    vm_ = collision(Erms.data[index], p);
                    a = dt * vm_ / 2;
                    alpha_ev = (1 - a) / (1 + a);
                    temp = 0.5 * e * dt * (1 + alpha_ev) / eps_0;
                } else if (niutype == 5) {
                    a = dt * Nu_c.data[index] / 2;
                    alpha_ev = (1 - a) / (1 + a);
                    temp = 0.5 * e * dt * (1 + alpha_ev) / eps_0;
                }
                Cevx.data[i * Cevx.ny + j] = temp * ne.data[index] / (1 + beta.data[index]);
            }
        }
        for (i = 0; i < Cevy.ny; ++i) {
            for (j = 0; j < Cevy.ny; j++) {
                index = m * i * beta.ny + m * j + m / 2;
                if (niutype == 4) {
                    vm_ = collision(Erms.data[index], p);
                    a = dt * vm_ / 2;
                    alpha_ev = (1 - a) / (1 + a);
                    temp = 0.5 * e * dt * (1 + alpha_ev) / eps_0;
                } else if (niutype == 5) {
                    a = dt * Nu_c.data[index] / 2;
                    alpha_ev = (1 - a) / (1 + a);
                    temp = 0.5 * e * dt * (1 + alpha_ev) / eps_0;
                }
                Cevy.data[i * Cevy.ny + j] = temp * ne.data[index] / (1 + beta.data[index]);
            }
        }
    }
    if (IsTEx) {
        for (i = 0; i < Cevz.nx; ++i)
            for (j = 0; j < Cevz.ny; j++) {
                index = m * i * beta.ny + m*j;
                if (niutype == 4) {
                    vm_ = collision(Erms.data[index], p);
                    a = dt * vm_ / 2;
                    alpha_ev = (1 - a) / (1 + a);
                    temp = 0.5 * e * dt * (1 + alpha_ev) / eps_0;
                } else if (niutype == 5) {
                    a = dt * Nu_c.data[index] / 2;
                    alpha_ev = (1 - a) / (1 + a);
                    temp = 0.5 * e * dt * (1 + alpha_ev) / eps_0;
                }
                Cevz.data[i * Cevz.ny + j] = temp * ne.data[index] / (1 + beta.data[index]);
            }
        //PrintData(Cevz);
    }
}

void InitVelocityCoeff_Equatioan_Five() {
    int i, j, index = 0;
    MyDataF temp, a = 0, gamma = 1;
    temp = e * dt * 0.5 / (me);
    if (IsTMx) {
        for (i = 0; i < Cvvx.nx; ++i) {
            for (j = 0; j < Cvvx.ny; j++) {
                index = (m * i + m / 2) * beta.ny + 2 * j;
                if (niutype == 4) {
                    a = dt * collision(Erms.data[index], p) / 2;
                    gamma = 1 + a;
                } else if (niutype == 5) {
                    a = dt * Nu_c.data[index] / 2;
                    gamma = 1 + a;
                }
                Cvvx.data[i * Cvvx.ny + j] = (1 - a) / gamma;
                Cvex.data[i * Cvex.ny + j] = temp / gamma;
            }
        }
        for (i = 0; i < Cvvy.ny; ++i) {
            for (j = 0; j < Cvvy.ny; j++) {
                index = m * i * beta.ny + m * j + m / 2;
                if (niutype == 4) {
                    a = dt * collision(Erms.data[index], p) / 2;
                    gamma = 1 + a;
                } else if (niutype == 5) {
                    a = dt * Nu_c.data[index] / 2;
                    gamma = 1 + a;
                }
                Cvvy.data[i * Cvvy.ny + j] = (1 - a) / gamma;
                Cvey.data[i * Cvey.ny + j] = temp / gamma;
            }
        }
    }
    if (IsTEx) {
        for (i = 0; i < Cvvz.nx; ++i) {
            for (j = 0; j < Cvvz.ny; j++) {
                index = m * i * beta.ny + m*j;
                if (niutype == 4) {
                    a = dt * collision(Erms.data[index], p) / 2;
                    gamma = 1 + a;
                } else if (niutype == 5) {
                    a = dt * Nu_c.data[index] / 2;
                    gamma = 1 + a;
                }
                Cvvz.data[i * Cvvz.ny + j] = (1 - a) / gamma;
                Cvez.data[i * Cvez.ny + j] = temp / gamma;
            }
        }
    }
}
/****************************************************************/
//coefficients when without density

void InitCeeNoDen() {
    int i, j;
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
    int i, j;
    //int index = 0;
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
        for (i = 0; i < Cehz.nx; ++i){
            for (j = 0; j < Cehz.ny; j++) {
                Cehz.data[i * Cehz.ny + j] = dt / (eps_0) / dx;
            }
        }
    }
    //PrintData(Cehz);
}

void InitCevNoDen() {
    int i, j;
    //int index = 0;
    //MyDataF temp;
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

/******************************************************************/
void UpdateCoeff() {

    if (IfWithDensity) {
        MyDataF gamma, a = 0;
        a = 0.5 * dt*vm;
        gamma = 1 + a;
        alpha = (1 - a) / (1 + a);

        InitBeta(&beta, gamma);
        /*PrintData(beta);*/
        InitCee(beta);
        InitCeh(beta);
        InitCev(beta, alpha);
        if (niutype == 4 || niutype == 5) {
            InitVelocityCoeff_Equatioan_Five();
        }
    }
};

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
        UpdateCoeff();
    } else {
        InitCeeNoDen();
        InitCehNoDen();
        InitCevNoDen();
    }
}
