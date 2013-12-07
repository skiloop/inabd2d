
/*
 *      This file defines functions that update the density.
 *
 *
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "fdtd.h"
#include "density.h"
#include "commonData.h"
#include "breakdownFormula.h"

void my_pause() {
    printf("Press any key to continue...");
    getchar();
}

void UpdateDensity() {
    int i, j, mt = m*tpis;
    MyDataF alpha_t, Eeff = 0, tau_m, kasi, Te;
    MyDataF ne_ij, neip1, neim1, nejm1, nejp1;
    //MyDataF temp2;
    MyDataF Deff;
    //MyDataF max_dtmui = 0, maxne = 0;
    //MyDataF mvi, mva, meff;
    //MyDataF nu_c;
    //MyDataF tmp = 0;
    MyDataF opt1, opt2, opt3;
    //MyDataF mopt1, mopt2, mopt3;

    //int ci = 0, cj = 0, mi = 0, mj = 0;
    int ind;
    MyDataF average_vi = 0, average_va = 0, max_vi = 0, max_va = 0;
    MyDataF dtF_Div_dsF2 = dt_F / ds_F / ds_F;
    int GridNum = (ne.nx - mt - mt)*(ne.ny - mt - mt);
    FILE *fp;
    MyDataF CErmsToEeff = 1.0 * pow(1.0 / (1.0 + omega * omega / vm / vm), 0.5) / 100.0;
    //PrintData(Ez_i);

    if (save_vi) {
        fp = fopen("vi.dat", "w");
        if (fp == NULL) {
            printf("Cannot open file vi.dat!\n");
            exit(-1);
        }
    }
    BackupMyStruct(ne, ne_pre);
#pragma omp parallel for num_threads(thread_count) schedule(dynamic) \
    private(i,j,k,ind, ne_ij, neip1, neim1, nejm1, nejp1,vi,va,opt1, opt2, opt3, \
    alpha_t, Eeff, tau_m, kasi, Te) //shared(Hx,Ez,Ey,pml,DA,DB,dy)
    for (i = mt; i < ne.nx - mt; i++) {
        for (j = mt; j < ne.ny - mt; j++) {
            ind = i * ne.ny + j;
            if (niutype != 4 && niutype != 5) {
                Eeff = Erms.data[ind] * CErmsToEeff;
            } else if (niutype == 5) {
                MyDataF tmp;
                tmp = Nu_c.data[ind] * Nu_c.data[ind];
                Eeff = Erms.data[ind] * sqrt(tmp / (tmp + omega * omega));
            }

            ne_ij = ne_pre.data[ind];
            neip1 = ne_pre.data[ind + ne.ny];
            neim1 = ne_pre.data[ind - ne.ny];
            nejm1 = ne_pre.data[ind - 1];
            nejp1 = ne_pre.data[ind + 1];
            if ((ne_ij <= 0 || fabs(Eeff) < 1e-10) && niutype != 4) {
                va = vi = 0;
            } else {
                switch (niutype) {
                    case 1:
                        Niu_MorrowAndLowke(&vi, &va, Eeff, NeutralGasDensityCM);
                        break;
                    case 2:
                        Niu_Nikonov(&vi, &va, Eeff, p);
                        break;
                    case 3:
                        Niu_Kang(&vi, &va, Eeff);
                        break;
                    case 4:
                        Te = ElectronTemperature(Erms.data[ind], p);
                        vi = (ne.data[ind] <= 0.0) ? 0 : ionization(Erms.data[ind], p);
                        Nu_c.data[ind] = collision(Erms.data[ind], p);
                        va = 0.0;
                        mu_e = e / me / Nu_c.data[ind]; //3.7e-2;
                        mu_i = mu_e / 100.0; //mu_e/mu_i ranges from 100 to 200//
                        De = mu_e*Te; //*1.6021e-19/e;//
                        Da = De * mu_i / mu_e;
                        break;
                    case 5:
                        Te = ElectronTemperature(Eeff, p);
                        if (ne.data[ind] <= 0.0) {
                            vi = 0.0;
                            va = 0.0;
                            rei = 0.0;
                        } else {
                            va = attachment(Eeff, p);
                            vi = ionization(Eeff, p);
                            rei = 1.0e-14;
                        }
                        Nu_c.data[ind] = collision(Eeff, p);
                        mu_e = e / me / Nu_c.data[ind]; //原来的是mu_e = e / me / vm; //
                        mu_i = mu_e / 100.0; //mu_e/mu_i ranges from 100 to 200//
                        De = mu_e*Te; //*1.6021e-19/e;//
                        Da = De * mu_i / mu_e;
                        break;
                    default:
                        alpha_t = Eeff / p;
                        if (alpha_t < 30) {
                            if (alpha_t < 1e-12) {
                                vi = 0;
                            } else if (alpha_t >= 1) {
                                vi = (1.45 + 0.01 * pow(alpha_t, 1.5))*2.5e7 * exp(-208 / alpha_t) * p;
                            } else {
                                vi = 5.14e11 * exp(-73 * pow(alpha_t, -0.44)) * p;
                            }
                        } else if (alpha_t > 120) {
                            if (alpha_t <= 3000) {
                                vi = 54.08e6 * pow(alpha_t, 0.5) * exp(-359 / alpha_t) * p;
                            } else {
                                vi = 5.14e11 * exp(-73 * pow(alpha_t, -0.44)) * p;
                            }
                        } else if (alpha_t > 54) {
                            vi = (1.32 + 0.054 * alpha_t)*1e7 * exp(-208 / alpha_t) * p;
                        } else {
                            vi = (5.0 + 0.19 * alpha_t)*1e7 * exp(-273.8 / alpha_t) * p;
                        }
                        //if(save_vi!=0)fprintf(fp,"%12f\t",vi);
                        va = 7.6e-4 * pow(alpha_t / (alpha_t + 218), 2) / p;
                }
            }

            if (ne_ij < 1) {
                Deff = De;
            } else {
                tau_m = eps_0 / (e * ne_ij * (mu_e + mu_i));
                kasi = vi * tau_m;
                Deff = (kasi * De + Da) / (kasi + 1);
            }

            opt1 = 1 + dt_F*vi;
            opt2 = Deff * (neip1 + neim1 + nejm1 + nejp1 - 4 * ne_ij) * dtF_Div_dsF2;
            opt3 = 1 + dt_F * (va + rei * ne_ij);
            ne.data[ind] = (ne_ij * opt1 + opt2) / opt3;

            // calculate average vi and va
            average_va += va / GridNum;
            average_vi += vi / GridNum;
#pragma omp critical
            {
                if (vi > max_vi)max_vi = vi;
            }
#pragma omp critical
            {
                if (va > max_va)max_va = va;
            }
        }
    }
    DensityBound(ne, (tpis-pis+1)*m, pis*m);
    printf("%5.4e\t%5.4e\t%5.4e\t%5.4e\t", average_va, average_vi, max_va, max_vi);
}

void InterpolatErms() {
    int l, k;
    int nk, nl;
    int ind;
    int is, js;
    int numx = 0, numy = 0;
    if (IsTEx) {
        numx = Ez_s.nx;
        numy = Ez_s.ny;
    } else if (IsTMx) {
        numx = Ey_s.nx;
        numy = Ex_s.ny;
    }
    is = tpis - 1;
    js = tpjs - 1;
    //for(k=0;k<Erms.nx-m;k++)
    //	for(l=0;l<Erms.ny-m;l++)
    //	{
    //		nk = k%m;
    //		nl = l%m;
    //		i = k/m;
    //		j = l/m;
    //		ind = i*Ez_s.ny+j;
    //		Erms.data[k*Erms.ny+l]=fabs((((m-nl)*(m-nk)*fabs(Ez_s.data[ind])+(m-nl)*nk*fabs(Ez_s.data[ind+Ez_s.ny])
    //		+nl*(m-nk)*fabs(Ez_s.data[ind+1])+(m-nl)*(m-nk)*fabs(Ez_s.data[ind+Ez_s.ny+1]))/m/m));
    //	//	if(k==cnepx && l== cnepy)
    //	//		printf("%4.4e\t%4.4e\t%4.4e\t%4.4e\n",Ez_s.data[ind],Ez_s.data[ind+1],Ez_s.data[ind+Ez_s.ny],Ez_s.data[ind+Ez_s.ny]);
    //	}
    for (k = is; k < numx - 1; k++) {
        for (l = js; l < numy - 1; l++) {
            ind = k * m * Erms.ny + l*m;
            for (nk = 1; nk <= m; nk++) {
                for (nl = 1; nl < m; nl++) {
                    Erms.data[ind + nk * Erms.ny + nl] = ((m - nl)*(m - nk) * Erms.data[ind]+(m - nl) * nk * Erms.data[ind + m * Erms.ny] +
                            nl * (m - nk) * Erms.data[ind + m] + nl * nk * Erms.data[ind + m * Erms.ny + m]) / m / m;
#ifdef _DEBUG
                    if (isnan(Erms.data[ind])) {
                        printf("Erms is nan at (%d,%d) at\n Line %d in function InterpolatErms()\n", nk*m, nl*m, __LINE__);
                    }
#endif
                }
            }
            for (nk = 1; nk < m; nk++) {
                Erms.data[ind + nk * Erms.ny + m] = ((m - nk) * Erms.data[ind + m] + nk * Erms.data[ind + m * Erms.ny + m]) / m;
#ifdef _DEBUG
                if (isnan(Erms.data[ind])) {
                    printf("Erms is nan at (%d,%d) at\n Line %d in function InterpolatErms()\n", l*m, nk*m, __LINE__);
                }
#endif
            }
        }
    }
    //printf("%4.4e\t%4.4e\t%4.4e\t%4.4e\n",Erms.data[cnepx*Erms.ny+cnepy],Erms.data[cnepx*Erms.ny+cnepy+m],
    //	Erms.data[cnepx*Erms.ny+cnepy+m*Erms.ny],Erms.data[cnepx*Erms.ny+cnepy+m*Erms.ny+m]);
}

void DensityBound(MyStruct stru, int bndwidth, const int swidth) {

    int i, j;
    int innerDown = swidth + bndwidth;
    int innerUp = stru.nx - 1 - innerDown;

    while (innerDown >= swidth) {
        int innerUpM1 = innerUp - 1;
        int innerUpM2 = innerUp - 2;
        int innerDownP1 = innerDown + 1;
        int innerDownP2 = innerDown + 2;

        // vertical lines
        for (j = innerDownP1; j <= innerUpM1; j++) {
            // left side line
            stru.data[innerDown * stru.ny + j] = 2 * stru.data[innerDownP1 * stru.ny + j] - stru.data[innerDownP2 * stru.ny + j];
            // right side line
            stru.data[innerUp * stru.ny + j] = 2 * stru.data[innerUpM1 * stru.ny + j] - stru.data[innerUpM2 * stru.ny + j];
        }

        // horizontal lines
        for (i = innerDownP1; i <= innerUpM1; i++) {
            // down side line
            stru.data[i * stru.ny + innerDown] = 2 * stru.data[i * stru.ny + innerDownP1] - stru.data[i * stru.ny + innerDownP2];
            // up side line
            stru.data[i * stru.ny + innerUp] = 2 * stru.data[i * stru.ny + innerUpM1] - stru.data[i * stru.ny + innerUpM2];
        }
        // Conner s
        //        stru.data[innerDown*stru.ny+innerDown]=2*stru.data[(innerDown+1)*stru.ny+innerDown+1]
        //                -stru.data[(innerDown+2)*stru.ny+innerDown+2];
        //        stru.data[innerDown*stru.ny+innerUp]=2*stru.data[(innerDown+1)*stru.ny+innerUp-1]
        //                -stru.data[(innerDown+2)*stru.ny+innerUp-2];
        //        stru.data[innerUp*stru.ny+innerDown]=2*stru.data[(innerUp-1)*stru.ny+innerDown+1]
        //                -stru.data[(innerUp-2)*stru.ny+innerDown+2];
        //        stru.data[innerUp*stru.ny+innerUp]=2*stru.data[(innerUp-1)*stru.ny+innerUp-1]
        //                -stru.data[(innerUp-2)*stru.ny+innerUp-2];
        // Conner s
        stru.data[innerDown * stru.ny + innerDown] = 0.5 * (
                2 * stru.data[(innerDownP1) * stru.ny + innerDown] - stru.data[(innerDownP2) * stru.ny + innerDown]
                + 2 * stru.data[(innerDown) * stru.ny + innerDownP1] - stru.data[(innerDown) * stru.ny + innerDownP2]);
        stru.data[innerDown * stru.ny + innerUp] = 0.5 * (
                2 * stru.data[(innerDownP1) * stru.ny + innerUp] - stru.data[(innerDownP2) * stru.ny + innerUp]
                + 2 * stru.data[(innerDown) * stru.ny + innerUpM1] - stru.data[(innerDown) * stru.ny + innerUpM2]);
        stru.data[innerUp * stru.ny + innerDown] = 0.5 * (
                2 * stru.data[(innerUpM1) * stru.ny + innerDown] - stru.data[(innerUpM2) * stru.ny + innerDown]
                + 2 * stru.data[(innerUp) * stru.ny + innerDownP1] - stru.data[(innerUp) * stru.ny + innerDownP2]);
        stru.data[innerUp * stru.ny + innerUp] = 0.5 * (
                2 * stru.data[(innerUp) * stru.ny + innerUpM1] - stru.data[(innerUp) * stru.ny + innerUpM2]
                + 2 * stru.data[(innerUpM1) * stru.ny + innerUp] - stru.data[(innerUpM2) * stru.ny + innerUp]);
        // increasment
        innerDown--;
        innerUp++;
    }
}

void CalErmsAtCoarseGrid() {
    int i, j, index;

    for (i = tpis; i <= tpie; i++) {
        for (j = tpjs; j <= tpje; j++) {
            index = i * m * Erms.ny + j*m;
            Erms.data[index] = sqrt(Erms.data[index] / MultiSize);
        }
    }
}

void CalSumESqrt() {
    int i, j, index, ind;
    if (IsTEx) {
        for (i = tpis; i <= tpie; i++) {
            for (j = tpjs; j <= tpje; j++) {
                index = i * Ez_s.ny + j;
                Erms.data[i * m * Erms.ny + j * m] += Ez_s.data[index] * Ez_s.data[index];
                //if(i == 141&&j==39)
                //	printf("%5.4e\t%5.4e\t%5.4e\t%5.4e\n",Erms.data[i*m*Erms.ny+j*m],Ez_s.data[index],Ez_s.data[index+Ez_s.ny]);
#ifdef _DEBUG
                if (isnan(Erms.data[index])) {
                    printf("Erms is nan at (%d,%d) at\n Line %d in function CalSumESqrt()\n", i*m, j*m, __LINE__);
                }
#endif
            }
        }
    } else if (IsTMx) {
        for (i = tpis; i <= tpie; i++) {
            for (j = tpjs; j <= tpje; j++) {
                ind = i * m * Erms.ny + j*m;
                index = i * Ex_s.ny + j;
                Erms.data[ind] += pow((Ex_s.data[index - Ex_s.nx] + Ex_s.data[index]) / 2, 2);
                index = i * Ey_s.ny + j;
                Erms.data[ind] += pow((Ey_s.data[index - 1] + Ey_s.data[index]) / 2, 2);
#ifdef _DEBUG
                if (isnan(Erms.data[ind])) {
                    printf("Erms is nan at (%d,%d) at\n Line %d in function CalSumESqrt()\n", i*m, j*m, __LINE__);
                }
#endif
            }
        }
    }
}

//算电场振幅

void CalSumESqrt_Emax() {
    int i, j, index, ind;
    MyDataF tmp = 0;
    if (IsTEx) {
        for (i = tpis; i <= tpie; i++) {
            for (j = tpjs; j <= tpje; j++) {
                index = i * Ez_s.ny + j;
                ind = (i * Erms.ny + j) * m;
                tmp = abs(Ez_s.data[index]);
                if (tmp > Erms.data[ind])
                    Erms.data[ind] = tmp;
            }
        }
    } else if (IsTMx) {
        for (i = tpis; i <= tpie; i++) {
            for (j = tpjs; j <= tpje; j++) {
                ind = i * m * Erms.ny + j*m;
                index = i * Ex_s.ny + j;
                tmp = pow((Ex_s.data[index - Ex_s.nx] + Ex_s.data[index]) / 2, 2);
                index = i * Ey_s.ny + j;
                tmp += pow((Ey_s.data[index - 1] + Ey_s.data[index]) / 2, 2);
                tmp = sqrt(abs(tmp));
                if (tmp > Erms.data[ind])
                    Erms.data[ind] = tmp;
                //if(i == 141&&j==39)
                //	printf("%5.4e\t%5.4e\t%5.4e\t%5.4e\n",Erms.data[i*m*Erms.ny+j*m],Ez_s.data[index],Ez_s.data[index+Ez_s.ny]);
            }
        }
    }
}

void calErmsAtCoarsePoint() {
    int i, j, hm;
    //int  sizex, sizey;
    //int index;
    int ind;
    MyDataF tmp;
    hm = m / 2;
    for (i = tpis; i <= tpie; i++) {
        for (j = tpjs; j <= tpje; j++) {
            ind = (i * m + hm) * Erms.ny + j*m;
            tmp = 0.25 * (Ey_s.data[i * Ey_s.ny + j] + Ey_s.data[(i - 1) * Ey_s.ny + j] + Ey_s.data[i * Ey_s.ny + j + 1] + Ey_s.data[(i - 1) * Ey_s.ny + j + 1]);
            Erms.data[ind] += Ex_s.data[i * Ex_s.ny + j] * Ex_s.data[i * Ex_s.ny + j] + tmp*tmp;
            ind = (i * m) * Erms.ny + j * m + hm;
            tmp = 0.25 * (Ex_s.data[i * Ex_s.ny + j] + Ex_s.data[(i + 1) * Ex_s.ny + j] + Ex_s.data[i * Ex_s.ny + j - 1] + Ex_s.data[(i + 1) * Ex_s.ny + j - 1]);
            Erms.data[ind] += Ey_s.data[i * Ey_s.ny + j] * Ey_s.data[i * Ey_s.ny + j] + tmp*tmp;
        }
    }
}

void calErmsAtCoarsePoint_Max() {
    int i, j;
    //int sizex, sizey;
    int hm;
    //int index;
    int ind;
    MyDataF tmp, Eabs = 0;
    hm = m / 2;
    for (i = tpis; i <= tpie; i++) {
        for (j = tpjs; j <= tpje; j++) {
            ind = (i * m + hm) * Erms.ny + j*m;
            tmp = 0.25 * (Ey_s.data[i * Ey_s.ny + j] + Ey_s.data[(i - 1) * Ey_s.ny + j] + Ey_s.data[i * Ey_s.ny + j + 1] + Ey_s.data[(i - 1) * Ey_s.ny + j + 1]);
            Erms.data[ind] += Ex_s.data[i * Ex_s.ny + j] * Ex_s.data[i * Ex_s.ny + j] + tmp*tmp;
            ind = (i * m) * Erms.ny + j * m + hm;
            tmp = 0.25 * (Ex_s.data[i * Ex_s.ny + j] + Ex_s.data[(i + 1) * Ex_s.ny + j] + Ex_s.data[i * Ex_s.ny + j - 1] + Ex_s.data[(i + 1) * Ex_s.ny + j - 1]);
            Eabs = abs(sqrt(Ey_s.data[i * Ey_s.ny + j] * Ey_s.data[i * Ey_s.ny + j] + tmp * tmp));
            if (Eabs > Erms.data[ind])Erms.data[ind] = Eabs;
        }
    }
}

void UpdateVelocity() {
    int i, j, index;
    //for format five and for
    if (niutype == 5 || niutype == 4) {
        if (IsTMx) {
            for (i = 0; i < Vex.nx; i++) {
                for (j = 0; j < Vex.ny; j++) {
                    index = i * Vex.ny + j;
                    Vex.data[index] = Cvvx.data[index] * Vex.data[index] - Cvex.data[index]*(Ex_s.data[index] + Ex_s_pre.data[index]);
                }
            }
            for (i = 0; i < Vey.nx; i++) {
                for (j = 0; j < Vey.ny; j++) {
                    index = i * Vey.ny + j;
                    Vey.data[index] = Cvvy.data[index] * Vey.data[index] - Cvey.data[index]*(Ey_s.data[index] + Ey_s_pre.data[index]);
                }
            }
        }
        if (IsTEx) {
            for (i = 0; i < Vez.nx; i++) {
                for (j = 0; j < Vez.ny; j++) {
                    index = i * Vez.ny + j;
                    Vez.data[index] = Cvvz.data[index] * Vez.data[index] - Cvez.data[index]*(Ez_s.data[index] + Ez_s_pre.data[index]);
                }
            }
        }

    } else {//for other format
        if (IsTMx) {
            for (i = 0; i < Vex.nx; i++) {
                for (j = 0; j < Vex.ny; j++) {
                    index = i * Vex.ny + j;
                    Vex.data[index] = alpha * Vex.data[index] - Cve * (Ex_s.data[index] + Ex_s_pre.data[ index]);
                }
            }
            for (i = 0; i < Vey.nx; i++) {
                for (j = 0; j < Vey.ny; j++) {
                    index = i * Vey.ny + j;
                    Vey.data[index] = alpha * Vey.data[index] - Cve * (Ey_s.data[ index] + Ey_s_pre.data[ index]);
                }
            }
        }
        if (IsTEx) {
            for (i = 0; i < Vez.nx; i++) {
                for (j = 0; j < Vez.ny; j++) {
                    index = i * Vez.ny + j;
                    Vez.data[index] = alpha * Vez.data[index] - Cve * (Ez_s_pre.data[index] + Ez_s.data[index]);
                }
            }
        }
    }
}
