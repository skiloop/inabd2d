
#include <math.h>

#include "common.h"
#include "dataType.h"
#include "commonData.h"
#include "fdtd.h"
#include "breakdownFormula.h"


int first, second, third;

void initCommonData() {
    /////////////////////////////////////
    //SET COMMOM DATA ZERO
    //////////////////////////////////////
    De = 0;
    Da = 0;
    //mu_e=0;
    mu_i = 0;
    vi = va = 0;
    kasi = 0;
    alpha = 0;
    //rei=0;


    //FDTD DATA
    dt = dt_F = dt_M = 0;
    ds_F = ds_M = 0;
    dx = 0;
    dy = 0;
    n = 0;

    //DOMAIN DATA
    nx = nxp1 = nxm1 = 0;
    ny = nyp1 = nym1 = 0;
    nbound = 0;

    //PARAMETERS OF INCIDENT WAVE
    f = 0; //frequency
    k = 0;
    E0 = H0 = 0;
    Hx0 = Hz0 = Hy0 = Ez0 = Ex0 = Ey0 = 0;
    lamda = 0;
    omega = 0;
    //////////////////////////////////////////////
    ///END OF SET COMMON DATA ZERO
    //////////////////////////////////////////////

    if (mu_e <= 0 || mu_e > 0.1)mu_e = e / me / vm; //3.7e-2;
    mu_i = mu_e / 100.0; //mu_e/mu_i ranges from 100 to 200//
    De = mu_e * 2 * 1.6021e-19 / e; //
    Da = De * mu_i / mu_e;
    vi = 5.4e9;
    va = 1.6e9;
    kasi = 0;
    if (rei < 0)rei = 0.0;
    D_kasi_max = (Da > De ? Da : De);

    f = FREQUENCY;
    lamda = c / f;
    E0 = E_0;
    H0 = E0 * sqrt(eps_0 / mu_0);
    k = omega / c;
    T = 1 / f;
    nbound = NUMBER_OF_CELLS_IN_PML_BOUND;
    omega = 2 * M_PI*f;

    Ratio_x = sin(phi); //sin(0.5*M_PI);
    Ratio_y = -cos(phi); //cos(0.5*M_PI);
    if (IsTMx) {
        Ex0 = Ratio_x*E0;
        Hz0 = H0;
        Ey0 = Ratio_y*E0;
    }
    if (IsTEx) {
        Hx0 = Ratio_x*H0;
        Hy0 = Ratio_y*H0;
        Ez0 = E0;
    }

    dx = dy = lamda / maxwellGridSize; //
    dt_M = dt = 0.5 * dx / c; //CourantFactor/c/sqrt(1/dx/dx+1/dy/dy+2.8e-13*M_PI*1e9);//
    //if(dt>1.0/12.0/f)dt=1.0/12.0/f;
    TotalTimeStep = (int) (0.5 + totaltime * 1e-9 / dt);
    SaveTimeStep = (int) (0.5 + 1e-9 / dt);
    if (SaveTimeStep < 0 && TotalTimeStep / SaveTimeStep > 10000) {
        SaveTimeStep = SAVE_LEAP;
    }
    ds_M = dx;
    ds_F = ds_M / m;
    dt_F = T; //CFL_factor*ds_F*ds_F*0.5/D_kasi_max;
    half_dt = dt / 2;


    /*if(dt_F*1000<dt_M||dt_F>dt_M)
            dt_F = dt_M/m;*/
    printf("\nT\t=\t%5.4e ns\nf\t=\t%-11lf GHz\n", T / 1e-9, f / 1e9);
    printf("mu_e\t=\t%5.4e\nmu_i\t=\t%5.4e\nDe\t=\t%5.4e\nDa\t=\t%5.4e\nrei\t=\t%5.4e\n", mu_e, mu_i, De, Da, rei);
    printf("lamda\t=\t%5.4e mm\nomega\t=\t%5.4e rad/s\n", lamda / 1e-3, omega);
    printf("dx\t=\t%5.4e mm\ndt\t=\t%5.4e ns\n", dx / 1e-3, dt / 1e-9);
    printf("\ndt_M\t=\t%5.4e ns\ndt_F\t=\t%5.4e ns\n", dt_M / 1e-9, dt_F / 1e-9);
    printf("ds_M\t=\t%5.4e mm\nds_F\t=\t%5.4e mm\n", ds_M / 1e-3, ds_F / 1e-3);
    printf("Ratio_x\t=\t%5.4e\nRatio_y\t=\t%5.4e\n", Ratio_x, Ratio_y);
    printf("\ndt_F/dt\t\t:\t%d\n", (int) (0.5 + dt_F / dt));
    printf("PML size\t:\t%d\n", nbound);
    printf("E0\t:\t%f\n", E0);
    printf("NE0\t:\t%f\n", NE0);
    printf("Fine Mesh Grid\t:\t%d\n", m);
    printf("Maxwell grid size\t:\t%d cells per lambda\n", maxwellGridSize);
    printf("Total propagating time\t:\t%f ns\n", totaltime);
    printf("niu function type:%d\n", niutype);
    printf("Save Time Step:%d\n", SaveTimeStep);
    printf("Total/Save:%d\n", TotalTimeStep / SaveTimeStep);
    printf("\n*************************************************************\n");
    //system("pause");
}

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
        for (i = 0; i < Cehz.nx; ++i) {
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
