/*
 *      This file defines functions that update the density.
 *
 *
 *
 *
 */
#include "InonizationFormula.h"

#include <stdlib.h>

#ifdef printf
#undef printf
#endif

void my_pause() {
    printf("Press any key to continue...");
    getchar();
}
void InterpolatErms();

void DensityBound(MyStruct stru, int bndwidth, const int);

void UpdateDensity() {
    int i, j, mt = m*tpis;
    MyDataF alpha_t, Eeff, tau_m, kasi;
    MyDataF ne_ij, temp2, neip1, neim1, nejm1, nejp1;
    MyDataF Deff;
    MyDataF max_dtmui = 0, maxne = 0, mvi, mva, meff;
    MyDataF tmp = 0;
    MyDataF opt1, opt2, opt3, mopt1, mopt2, mopt3;
    int ci = 0, cj = 0, mi = 0, mj = 0;
    int ind;
    MyDataF average_vi = 0, average_va = 0, max_vi = 0, max_va = 0;
    int GridNum = (ne.nx - mt - mt)*(ne.ny - mt - mt);
    FILE *fp;
    //PrintData(Ez_i);

    if (save_vi != 0) {
        fp = fopen("vi.dat", "w");
        if (fp == NULL) {
            printf("Cannot open file vi.dat!\n");
            exit(-1);
        }
    }
    BackupMyStruct(ne, ne_pre);
    //printf("%5.4e\t%5.4e\t",Erms.data[250*Erms.ny+250],ne.data[250*ne.ny+250]);
    for (i = mt; i < ne.nx - mt; i++)
        for (j = mt; j < ne.ny - mt; j++) {
            ind = i * ne.ny + j;
            Eeff = Erms.data[ind] / 100 * pow(1 / (1 + omega * omega / vm / vm), 0.5);

            ne_ij = ne_pre.data[ind];
            neip1 = ne_pre.data[ind + ne.ny];
            neim1 = ne_pre.data[ind - ne.ny];
            nejm1 = ne_pre.data[ind - 1];
            nejp1 = ne_pre.data[ind + 1];
            if (ne_ij < 1 || fabs(Eeff) < 0.1) va = vi = 0;
            else {
                switch (niutype) {
                    case 1:
                        Niu_MorrowAndLowke(&vi, &va, Eeff, ne_ij * 1e6);
                        break;
                    case 2:
                        Niu_Nikonov(&vi, &va, Eeff, p);
                        break;
                    case 3:
                        Niu_Kang(&vi, &va, Eeff);
                        break;
                    default:
                        alpha_t = Eeff / p;
                        if (alpha_t < 30) {
                            if (alpha_t < 1e-12) {
                                vi = 0;
                            } else {
                            }
                        } else if (alpha_t > 120) {
                            vi = 54.08e6 * pow(alpha_t, 0.5) * exp(-359 / alpha_t) * p;
                        } else if (alpha_t <= 54) {
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
            /*temp2=ne.data[i*ne.ny+j]=A*pow(CTime,-1.5)*exp(vi*CTime-(dx*dx+dy*dy)/(4*Deff*CTime));
             */


            opt1 = 1 + dt_F*vi;
            opt2 = Deff * dt_F * (neip1 + neim1 + nejm1 + nejp1 - 4 * ne_ij) / ds_F / ds_F;
            opt3 = 1 + dt_F * (va + rei * ne_ij);
            temp2 = ne.data[ind] = (ne_ij * opt1 + opt2) / opt3;

            // calculate average vi and va
            average_va += va / GridNum;
            average_vi += vi / GridNum;
            if (vi > max_vi)max_vi = vi;
            if (va > max_va)max_va = va;
            if (maxne < ne.data[ind]) {
                maxne = ne.data[ind];
                mi = i;
                mj = j;
                mopt1 = opt1;
                mopt2 = opt2;
                mopt3 = opt3;
                mvi = vi;
                mva = va;
                meff = Eeff;
                //#ifdef _DEBUG
                //				if (isnan(maxne) != 0)
                //				{
                //					printf("maxne is nan!");
                //					my_pause();
                //				}
                //				if (isnan(opt1) != 0)
                //				{
                //					printf("opt1 is nan!");
                //					getchar();
                //				}
                //				if (isnan(opt2) != 0)
                //				{
                //					printf("opt2 is nan!");
                //					my_pause();
                //				}
                //				if (isnan(opt3) != 0)
                //				{
                //					printf("opt4 is nan!");
                //					my_pause();
                //				}
                //				if (isnan(mvi) != 0)
                //				{
                //					printf("mvi is nan!");
                //					my_pause();
                //				}
                //				if (isnan(mva) != 0)
                //				{
                //					printf("mva is nan!");
                //					my_pause();
                //				}
                //				if (isnan(meff) != 0)
                //				{
                //					printf("meff is nan!");
                //					my_pause();
                //				}
                //#endif
            }

            /*
            temp2=ne.data[ind]=
                        (ne_ij*(1+dt_F*vi))
                        /(1+dt_F*(va+rei*ne_ij));
             */
            /*
            if(i==cnepx&&j==cnepy){
                //printf("%4.4e\n",Erms.data[ind]);

                printf("\n\nvi\t=\t%6.5e\ntau\t=\t%6.5e\n",vi,tau_m);
                printf("dt_f*vi\t=\t%6.5lf\nalpha\t=\t%5.5e\nEeff\t=\t%5.5e\nErms\t=\t%4.4e\n",dt_F*vi,alpha_t,Eeff,Erms.data[ind]);
                printf("Deff\t=\t%4.3e\nkasi\t=\t%4.3e\n",Deff,kasi);
                printf("De\t=\t%4.3e\nDa\t=\t%5.4e\n\n",De,Da);
                printf("Sp\t\t=\t%6.5e\n",Deff*dt_F*(neip1+neim1+nejm1+nejp1-4*ne_ij)/ds_F/ds_F);
                printf("Neij_dt_vi\t=\t%6.5e\n",ne_ij*dt_F*vi);
                printf("%5.4e\t%5.4e\t%5.4e\t%5.4e\n%5.4e\n",neip1,neim1,nejm1,nejp1,ne_ij);
                printf("(1+dt_F*(va+rei*ne_ij)\t=\t%5.4e\n",(1+dt_F*(va+rei*ne_ij)));
                printf("Deff*dtF/dsf2 =\t%6.5e\n",Deff*dt_F/ds_F/ds_F);
                printf("Nes\t=\t%6.5e\n",neip1+neim1+nejm1+nejp1-4*ne_ij);
                printf("neij\t=\t%6.5e\n",ne_ij);
                printf("ne4pm\t=\t%6.5e\n",neip1+neim1+nejm1+nejp1);

                printf("ne_new\t=\t%6.5e\n",ne.data[ind]);
                printf("increased\t=\t%6.5e\n",ne_ij*dt_F*vi+Deff*dt_F*(neip1+neim1+nejm1+nejp1-4*ne_ij)/ds_F/ds_F);
                printf("\n**************************************************\n");

            }
             */
            //tmp = dt_F*vi;//(1+dt_F*vi)/(1+dt_F*va);//Deff*dt_F/ds_F/ds_F/(1+dt_F*va);
            if (opt1 > max_dtmui) {
                //printf("max_dtmui = %f\n",max_dtmui);
                ci = i;
                cj = j;
                max_dtmui = opt1;
                //printf("Deff*dtf/ds/ds/(1+dtf*va) = %5.3e\n",Deff*dt_F/ds_F/ds_F/(1+dt_F*va));
                //puts("\n*******************************************");
                //printf("Eeff = %4.2e\nvi = %5.2e\nalphat = %4.2e\n",Eeff,vi,alpha_t);
                //printf("Step = %d\ndmax_dtmui = %5.2e\ndtf*vi = %4.2e\nva = %4.2e\n",Step,max_dtmui,dt_F*vi,va);
                //printf("Ez_s(%d,%d) = %3.3e\tne = %4.2e\n",i/m,j/m,Ez_s.data[(i/m)*Ez.ny+j/m],temp2);
                //CaptData(0,"Ezs.dat",Ez_s);
                //CaptData(0,"ne.dat",ne);
                ////puts("*******************************************");
                //system("cls");
            }
            /*
#ifdef DEBUG
            if((temp2>MAX_NE)||(temp2<-MAX_NE))
            {
                printf("\n==========================================================\n");
                printf("Deff\t=\t%5.2e\n",Deff);
                printf("Eeff\t=\t%5.2e\nvi\t=\t%5.2e\nalphat\t=\t%5.2e\n",Eeff,vi,alpha_t);
                printf("Step\t=\t%d\ndmax_dtmui\t=\t%5.2e\ndtf*vi\t=\t%5.2e\nva\t=\t%5.2e\n",Step,max_dtmui,dt_F*vi,va);
                printf("Erms(%d,%d)\t=%5.2e\tne\t=%5.2e\n",i,j,Erms.data[i*Erms.ny+j],temp2);
                printf("\n==========================================================\n");
                //	CaptData(0,"Ezs.dat",Ez_s);
                //	CaptData(0,"ne.dat",ne);
                //	//system("pause");
                ne.data[i*ne.ny+j] = MAX_NE;
            }
#endif
             */
            /*
#ifdef _DEBUG
            if (isnan(ne.data[ind]))
            {
                printf("ne is nan at (%d,%d) at\n Line %d in function UpdateDensity()\n",i,j,__LINE__);
            }
#endif
             */

        }

    //DensityBound(ne,m*tpis,0);
    //PrintData(ne);
    //if(max_dtmui>1.0)
    //    printf("%5.4e\t%5.4e\t",Erms.data[ci*Erms.ny+cj],ne.data[ci*ne.ny+cj]);
    //printf("%5f\t",max_dtmui);
    //printf("%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t",maxne,meff,mvi,mva,mopt1,mopt2,mopt3,Erms.data[mi*Erms.ny+mj]);
    //fprintf(filedat,"%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t%5.4e\t",maxne,meff,mvi,mva,mopt1,mopt2,mopt3,Erms.data[mi*Erms.ny+mj]);
    //fprintf(filedat,"%5.4e\t%5.4e\t",Erms.data[ci*Erms.ny+cj],ne.data[ci*ne.ny+cj]);
    //fprintf(filedat,"%5f\n",max_dtmui);
    //
    //
    //
    //system("pause");
    printf("%5.4e\t%5.4e\t%5.4e\t%5.4e\t", average_va, average_vi, max_va, max_vi);

}

void InterpolatErms() {
    int l, k;
    int nk, nl;
    int ind;
    int numx, numy;
    if (IsTEx) {
        numx = Ez_s.nx;
        numy = Ez_s.ny;
    } else if (IsTMx) {
        numx = Ey_s.nx;
        numy = Ex_s.ny;
    }
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
    for (k = 0; k < numx - 1; k++)
        for (l = 0; l < numy - 1; l++) {
            ind = k * m * Erms.ny + l*m;
            for (nk = 1; nk <= m; nk++)
                for (nl = 1; nl < m; nl++) {
                    Erms.data[ind + nk * Erms.ny + nl] = ((m - nl)*(m - nk) * Erms.data[ind]+(m - nl) * nk * Erms.data[ind + m * Erms.ny] +
                            nl * (m - nk) * Erms.data[ind + m] + nl * nk * Erms.data[ind + m * Erms.ny + m]) / m / m;
#ifdef _DEBUG
                    if (isnan(Erms.data[ind])) {
                        printf("Erms is nan at (%d,%d) at\n Line %d in function InterpolatErms()\n", nk*m, nl*m, __LINE__);
                    }
#endif
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
    //printf("%4.4e\t%4.4e\t%4.4e\t%4.4e\n",Erms.data[cnepx*Erms.ny+cnepy],Erms.data[cnepx*Erms.ny+cnepy+m],
    //	Erms.data[cnepx*Erms.ny+cnepy+m*Erms.ny],Erms.data[cnepx*Erms.ny+cnepy+m*Erms.ny+m]);
}

void DensityBound(MyStruct stru, int bndwidth, const int swidth) {

    int i, j, ind;
    int ii, i1, i2;
    double tmp;
    ii = swidth + bndwidth;
    i1 = stru.nx - swidth - bndwidth;
    i2 = stru.nx - swidth;
    //Section 1
    for (j = swidth; j < i1; j++) {
        ind = ii * stru.ny + j;
        tmp = 0; //2*stru.data[ind+stru.ny]-stru.data[ind+2*stru.ny];
        for (i = swidth; i <= ii; i++)
            stru.data[i * stru.ny + j] = tmp;
    }
    //Section 2
    for (i = 0; i < i1; i++) {
        ind = i * stru.ny + i1;
        tmp = 0; //2*stru.data[ind-1]-stru.data[ind-2];
        for (j = i2 - 1; j >= i1; j--)
            stru.data[i * stru.ny + j] = tmp;
    }
    //Section 3
    for (j = ii; j < i2; j++) {
        ind = i1 * stru.ny + j;
        tmp = 0; //2*stru.data[ind-stru.ny]-stru.data[ind-2*stru.ny];
        for (i = i2 - 1; i >= i1; i--)
            stru.data[i * stru.ny + j] = tmp;
    }

    //Section 4
    for (i = i1; i < i2; i++) {
        ind = ii + i * stru.ny;
        tmp = 0; //2*stru.data[ind+1]-stru.data[ind+2];
        for (j = swidth; j <= ii; j++)
            stru.data[i * stru.ny + j] = tmp;
    }
}

void CalErmsAtCoarseGrid() {
    int i, j, index;

    for (i = pis - 1; i <= pie + 1; i++)
        for (j = pjs - 1; j <= pje + 1; j++) {
            index = i * m * Erms.ny + j*m;
            Erms.data[index] = sqrt(Erms.data[index] / MultiSize);
            //if(Erms.data[i*m*Erms.ny+j*m]>E0)
            //{
            //	printf("i = %d\tj = %d\nErms = %6.4e\n",i,j,Erms.data[i*m*Erms.ny+j*m]);

            //}
            //if(i ==cnepx && j == cnepy)//(i <= cnepx +m)&&(i >= cnepx) && (j <= cnepy+m) &&(j >= cnepy))
            //		printf("Erms = %5f\n",Erms.data[i][j]);
#ifdef _DEBUG
            if (isnan(Erms.data[index])) {
                printf("Erms is nan at (%d,%d) at\n Line %d in function CalErmsAtCoarseGrid()\n", i*m, j*m, __LINE__);
            }
#endif

        }
}

//

void CalSumESqrt() {
    int i, j, index, ind;
    if (IsTEx) {
        for (i = tpis; i <= tpie; i++)
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
    } else if (IsTMx) {
        for (i = tpis; i <= tpie; i++)
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
                //if(i == 141&&j==39)
                //	printf("%5.4e\t%5.4e\t%5.4e\t%5.4e\n",Erms.data[i*m*Erms.ny+j*m],Ez_s.data[index],Ez_s.data[index+Ez_s.ny]);
            }
    }
}
