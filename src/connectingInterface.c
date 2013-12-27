
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include<engine.h>
#include<mex.h>
#ifdef printf
#undef printf
#endif
#include "common.h"
#include "fdtd.h"
#include "dataType.h"
#include "commonData.h"
#include "connectingInterface.h"


MyDataF *Ei, *Hi;
int eilen, hilen;
MyDataF CMur;

MyDataF Ceyhz;
MyDataF Cexhz;
MyDataF chiei, ceihi;
MyDataF Vp_ratio; //phase velocity ratio


//Hy Delay
MyDataF *Dhyl, *Dhyr;
//Hx Delay
MyDataF *Dhxb, *Dhxt;
//Ez Delay
MyDataF *Dezb, *Dezt, *Dezl, *Dezr;

//Ey Delay
MyDataF *Deyl, *Deyr;
//Ex Delay
MyDataF *Dexb, *Dext;
//Hz Delay
MyDataF *Dhzb, *Dhzt, *Dhzl, *Dhzr;

MyDataF t_per_cell; //time wave propagates per Yee cell
int xs, ys, xe, ye; //
MyDataF ds;

mxArray *myArr;
extern Engine* ep;
#define INC_SIZE 0

void initSource() {
    switch (srcType) {
        case GAUSSIAN:
            tauSquare = T * T / 4;
            t_0 = T * 1.2;
            pSource = GaussianSource;
            break;
        case COSINE:
            pSource = CosineSource;
            break;
        case SINE:
        default:
            pSource = SineSource;
    }
}

void initconnect() {
    //#ifdef  MATLAB_SIMULATION
    //	myArr *MyArray;
    //
    //#endif
    int i;
    MyDataF PhaseVelRatio(MyDataF);

    eilen = (int) ceil(sqrt((ye - ys)*(ye - ys)+(xe - xs)*(xe - xs))) + 5;
    hilen = eilen - 1;

    //ds = dx * PhaseVelRatio(0) / PhaseVelRatio(phi);
    ds=dx;
    CMur = (c * dt - ds) / (c * dt + ds);

    if (IsTMx) {
        chiei = dt / mu_0 / ds;
        ceihi = dt / eps_0 / ds;
    }
    if (IsTEx) {
        chiei = dt / mu_0 / ds;
        ceihi = dt / eps_0 / ds;
    }


    Ei = (MyDataF*) malloc(eilen * sizeof (MyDataF));
    Hi = (MyDataF*) malloc(hilen * sizeof (MyDataF));

    if (Ei == NULL || Hi == NULL) {
        exit(EXIT_FAILURE);
    }

    for (i = 1; i < eilen; i++)Ei[i] = 0;
    for (i = 0; i < hilen; i++)Hi[i] = 0;
    Ei[0] = E0 * Source(0);
//#ifdef  MATLAB_SIMULATION
//    myArr = mxCreateDoubleMatrix(1, eilen, mxREAL);
//    memcpy(mxGetPr(myArr), Ei, eilen * sizeof (MyDataF));
//    engPutVariable(ep, "Ei", myArr);
//    engEvalString(ep, "figure,eih=plot(Ei);grid on;");
//    mxDestroyArray(myArr);
//#endif
}

void econnect(MyDataF t) {
    //#ifdef  MATLAB_SIMULATION
    //	myArr *MyArray;
    //
    //#endif
    int i, index, ind;
    int di;

    MyDataF ei_last2, df;
    ei_last2 = Ei[eilen - 2];
    Ei[0] = E0 * Source(t);
    for (i = 1; i <= eilen - 2; i++) {
        Ei[i] += ceihi * (Hi[i] - Hi[i - 1]);
    }

    //mur boundary 
    Ei[eilen - 1] = ei_last2 + CMur * (Ei[eilen - 2] - Ei[eilen - 1]);

    //#ifdef  MATLAB_SIMULATION
    //    myArr = mxCreateDoubleMatrix(1, eilen, mxREAL);
    //    memcpy(mxGetPr(myArr), Ei, eilen * sizeof (MyDataF));
    //    engPutVariable(ep, "Ei", myArr);
    //    engEvalString(ep, "set(eih,'ydata',Ei);"); //
    //    //engEvalString(ep,"title(ind);ind=ind+1;");
    //    mxDestroyArray(myArr);
    //#endif
    if (IsTEx) {
        for (ind = tpis, i = 0; ind <= tpie; ind++, i++) {
            index = ind * Ez.ny;
            //bottom 
            di = (int) floor(Dhxb[i]);
            df = Dhxb[i] - di;
            Ez.data[index + tpjs] += Cehz.data[index + tpjs] * RatioHx * (Hi[di] + df * (Hi[di + 1] - Hi[di]));

            //top
            di = (int) floor(Dhxt[i]);
            df = Dhxt[i] - di;
            Ez.data[index + tpje] -= Cehz.data[index + tpje] * RatioHx * (Hi[di] + df * (Hi[di + 1] - Hi[di]));
        }
        //bound xn,xp
        for (ind = tpjs, i = 0; ind <= tpje; ind++, i++) {
            //left side
            index = tpis * Ez.ny + ind;
            di = (int) floor(Dhyl[i]);
            df = Dhyl[i] - di;
            Ez.data[index] -= Cehz.data[index] * RaitoHy * (Hi[di] + df * (Hi[di + 1] - Hi[di]));

            //right side
            index = tpie * Ez.ny + ind;
            di = (int) floor(Dhyr[i]);
            df = Dhyr[i] - di;
            Ez.data[index] += Cehz.data[index] * RaitoHy * (Hi[di] + df * (Hi[di + 1] - Hi[di]));
        }
    }
    if (IsTMx) {
        int i, j;
        int ind;
        for (ind = 0, j = tpjs; j < tpje; j++, ind++) {
            di = (int) floor(Dhzl[ind]);
            df = Dhzl[ind] - di;
            Ey.data[tpis * Ey.ny + j] -= Ceyhz * RatioHz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //left side

            di = (int) floor(Dhzr[ind]);
            df = Dhzr[ind] - di;
            Ey.data[tpie * Ey.ny + j] += Ceyhz * RatioHz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //right side
        }
        for (ind = 0, i = tpis; i < tpie; i++, ind++) {
            di = (int) floor(Dhzb[ind]);
            df = Dhzb[ind] - di;
            Ex.data[i * Ex.ny + tpjs] -= Cexhz * RatioHz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //bottom side

            di = (int) floor(Dhzt[ind]);
            df = Dhzt[ind] - di;
            Ex.data[i * Ex.ny + tpje] += Cexhz * RatioHz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //top side
        }
    }
}

void mconnect(MyDataF t) {

    int i, ind, di;
    MyDataF df;

    for (i = 0; i < hilen; i++) {
        Hi[i] = Hi[i] + chiei * (Ei[i + 1] - Ei[i]);
    }

    int mtpis, mtpjs;
    mtpis = tpis - 1;
    mtpjs = tpjs - 1;
    if (IsTMx) {
        for (ind = tpis, i = 0; ind < tpie; ind++, i++) {
            di = (unsigned) floor(Dexb[i]);
            df = Dexb[i] - di;
            Hz.data[ind * Hz.ny + mtpjs] -= Chzex * RatioEx * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //bottom side
            di = (unsigned) floor(Dext[i]);
            df = Dext[i] - di;
            Hz.data[ind * Hz.ny + tpje] += Chzex * RatioEx * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //top side
        }
        for (ind = tpjs, i = 0; ind < tpje; ind++, i++) {
            di = (unsigned) floor(Deyl[i]);
            df = Deyl[i] - di;
            Hz.data[mtpis * Hz.ny + ind] -= Chzey * RatioEy * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //left side
            di = (unsigned) floor(Deyr[i]);
            df = Deyr[i] - di;
            Hz.data[tpie * Hz.ny + ind] += Chzey * RatioEy * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //right side
        }
    }
    if (IsTEx) {
        for (ind = tpjs, i = 0; ind <= tpje; ind++, i++) {
            //left side
            di = (int) floor(Dezl[i]);
            df = Dezl[i] - di;
            Hy.data[mtpis * Hy.ny + ind] -= Chyez * RatioEz * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //0;//
            //right side
            di = (int) floor(Dezr[i]);
            df = Dezr[i] - di;
            Hy.data[tpie * Hy.ny + ind] += Chyez * RatioEz * (Ei[di] + df * (Ei[di + 1] - Ei[di]));
        }
        //Adjust Hy
        for (ind = tpis, i = 0; ind <= tpie; ind++, i++) {
            // bottom
            di = (int) floor(Dezb[i]);
            df = Dezb[i] - di;
            Hx.data[ind * Hx.ny + mtpjs] -= Chxez * RatioEz * (Ei[di] + df * (Ei[di + 1] - Ei[di]));
            //yp
            di = (int) floor(Dezt[i]);
            df = Dezt[i] - di;
            Hx.data[ind * Hx.ny + tpje] += Chxez * RatioEz * (Ei[di] + df * (Ei[di + 1] - Ei[di]));
        }
    }

}

void FreeDelayArrays() {
    if (IsTEx) {
        free(Dezb);
        free(Dezt);
        free(Dezl);
        free(Dezr);

        free(Dhxb);
        free(Dhxt);
        free(Dhyl);
        free(Dhyr);
    }
    if (IsTMx) {
        free(Dhzb);
        free(Dhzt);
        free(Dhzl);
        free(Dhzr);

        free(Dexb);
        free(Dext);
        free(Deyl);
        free(Deyr);
    }
}

void end_connect() {
    free(Ei);
    free(Hi);
    FreeDelayArrays();
}

void CalDelays() {
    int i;
    int Max_Index_Dez_x;
    int Max_Index_Dez_y;

    int Max_Index_Dhx;
    int Max_Index_Dhy;

    MyDataF delays(MyDataF px, MyDataF py);


    t_per_cell = dx / c;

    // normalize phi
    phi = phi - 2.0 * M_PI * (floor(phi / (2.0 * M_PI)));

    int mtpis = tpis - 1;
    int mtpjs = tpjs - 1;
    int mtpie = tpie + 1;
    int mtpje = tpje + 1;
    MyDataF mei = 1.0; // distance from Om to O0 in ds
    MyDataF mhi = mei - 0.5;
    // deside original point
    if ((phi <= 0.5 * M_PI)) {
        xs = mtpis;
        ys = mtpjs;
        xe = mtpie;
        ye = mtpje;
    } else if ((phi <= M_PI)) {
        xs = mtpie;
        ys = mtpjs;
        xe = mtpis;
        ye = mtpje;
    } else if ((phi <= 1.5 * M_PI)) {
        xs = mtpie;
        ys = mtpje;
        xe = mtpis;
        ye = mtpjs;
    } else {
        xs = mtpis;
        ys = mtpje;
        xe = mtpie;
        ye = mtpjs;
    }
    if (IsTEx) {
        Max_Index_Dhx = Max_Index_Dez_x = abs(xe - xs);
        Max_Index_Dhy = Max_Index_Dez_y = abs(ye - ys);

        Dhxb = (MyDataF*) malloc((Max_Index_Dhx + 1) * sizeof (MyDataF));
        Dhxt = (MyDataF*) malloc((Max_Index_Dhx + 1) * sizeof (MyDataF));
        if (Dhxb == NULL || Dhxt == NULL) {
            printf("Cannot allocate enough memory for Dhx!\n");
            exit(EXIT_FAILURE);
        }
        Dhyl = (MyDataF*) malloc((Max_Index_Dhy + 1) * sizeof (MyDataF));
        Dhyr = (MyDataF*) malloc((Max_Index_Dhy + 1) * sizeof (MyDataF));
        if (Dhyl == NULL || Dhyr == NULL) {
            printf("Cannot allocate enough memory for Dhy!\n");
            exit(EXIT_FAILURE);
        }

        Dezb = (MyDataF*) malloc((Max_Index_Dez_x + 1) * sizeof (MyDataF));
        Dezt = (MyDataF*) malloc((Max_Index_Dez_x + 1) * sizeof (MyDataF));
        Dezl = (MyDataF*) malloc((Max_Index_Dez_y + 1) * sizeof (MyDataF));
        Dezr = (MyDataF*) malloc((Max_Index_Dez_y + 1) * sizeof (MyDataF));

        if (Dezt == NULL || Dezb == NULL || Dezl == NULL || Dezr == NULL) {
            printf("Cannot allocate enough memory for Dez!\n");
            exit(EXIT_FAILURE);
        }

        for (i = 0; i <= Max_Index_Dhx; i++) {
            Dhxb[i] = delays(tpis + i, tpjs - 0.5) + mhi;
            Dhxt[i] = delays(tpis + i, tpje + 0.5) + mhi;
        }
        for (i = 0; i <= Max_Index_Dhy; i++) {
            Dhyl[i] = delays(tpis - 0.5, tpjs + i) + mhi;
            Dhyr[i] = delays(tpie + 0.5, tpjs + i) + mhi;
        }
        for (i = 0; i <= Max_Index_Dez_x; i++) {
            Dezb[i] = delays(tpis + i, tpjs) + mei;
            Dezt[i] = delays(tpis + i, tpje) + mei;
        }
        for (i = 0; i <= Max_Index_Dez_y; i++) {
            Dezl[i] = delays(tpis, tpjs + i) + mei;
            Dezr[i] = delays(tpie, tpjs + i) + mei;
        }
    }
    if (IsTMx) {
        //Ex Ey are total fields and Hz is scattered field
        int Max_Index_Dhz_x;
        int Max_Index_Dhz_y;

        int Max_Index_Dex;
        int Max_Index_Dey;

        Max_Index_Dex = Max_Index_Dhz_x = (int) abs((MyDataF) (xe - xs));
        Max_Index_Dey = Max_Index_Dhz_y = (int) abs((MyDataF) (ye - ys));

        Dexb = (MyDataF*) calloc((Max_Index_Dex), sizeof (MyDataF));
        Dext = (MyDataF*) calloc((Max_Index_Dex), sizeof (MyDataF));

        if (Dexb == NULL || Dext == NULL) {
            printf("Cannot allocate enough memory for Dex!");
            exit(EXIT_FAILURE);
        }
        Deyl = (MyDataF*) calloc((Max_Index_Dey), sizeof (MyDataF));
        Deyr = (MyDataF*) calloc((Max_Index_Dey), sizeof (MyDataF));
        if (Deyl == NULL || Deyr == NULL) {
            printf("Cannot allocate enough memory for Dey!\n");
            ;
            exit(EXIT_FAILURE);
        }

        Dhzb = (MyDataF*) calloc((Max_Index_Dhz_x), sizeof (MyDataF));
        Dhzt = (MyDataF*) calloc((Max_Index_Dhz_x), sizeof (MyDataF));
        Dhzl = (MyDataF*) calloc((Max_Index_Dhz_y), sizeof (MyDataF));
        Dhzr = (MyDataF*) calloc((Max_Index_Dhz_y), sizeof (MyDataF));
        if (Dhzt == NULL || Dhzb == NULL || Dhzl == NULL || Dhzr == NULL) {
            fprintf(stderr, "Cannot allocate enough memory for Dhz!");
            exit(EXIT_FAILURE);
        }

        for (i = 0; i < Max_Index_Dex; i++) {
            Dexb[i] = delays(tpis + i + 0.5, tpjs) + mei;
            Dext[i] = delays(tpis + i + 0.5, tpje) + mei;
        }
        for (i = 0; i < Max_Index_Dey; i++) {
            Deyl[i] = delays(tpis, tpjs + i + 0.5) + mei;
            Deyr[i] = delays(tpie, tpjs + i + 0.5) + mei;
        }
        for (i = 0; i < Max_Index_Dhz_x; i++) {
            Dhzb[i] = delays(tpis + i + 0.5, tpjs - 0.5) + mhi;
            Dhzt[i] = delays(tpis + i + 0.5, tpje + 0.5) + mhi;
        }
        for (i = 0; i < Max_Index_Dhz_y; i++) {
            Dhzl[i] = delays(tpis - 0.5, tpjs + i + 0.5) + mhi;
            Dhzr[i] = delays(tpie + 0.5, tpjs + i + 0.5) + mhi;
        }
    }

}

void IntTtlFldDmnBnd() {
    tpis = pis + SCATTER_FIELD_DOMAIN_BND_SIZE;
    tpjs = pjs + SCATTER_FIELD_DOMAIN_BND_SIZE;
    tpie = pie - SCATTER_FIELD_DOMAIN_BND_SIZE;
    tpje = pje - SCATTER_FIELD_DOMAIN_BND_SIZE;
    printf("\n*************************************************************\n");
    printf("\ntpis\t=\t%d\ntpjs\t=\t%d\ntpie\t=\t%d\ntpje\t=\t%d\n", tpis, tpjs, tpie, tpje);
    CalDelays();
}

MyDataF PhaseVelRatio(MyDataF angle) {
    MyDataF A, B, C, S, k, kp;
    const int maxIter = 100000;
    int iter = 0;
    S = c * dt / dx;
    A = dx * cos(angle);
    B = dx * sin(angle);
    C = pow(sin(M_PI * S / maxwellGridSize) / S, 2.0);
    C = 2 * (C - 1);
    k = 2 * M_PI;
    kp = 0.2;
    while (fabs(k - kp) > 1e-5 && iter < maxIter) {
        kp = k;
        k = k + (cos(A * k) + cos(B * k) + C) / (A * sin(A * k) + B * sin(B * k));
        iter++;
    }
    if (maxIter==iter){
        k=k>kp?kp:k;
    }
    return 2 * M_PI * c / k;
}

MyDataF GaussianSource(MyDataF t) {

    return exp(-(t - t_0)*(t - t_0) / tauSquare);
}

MyDataF SineSource(MyDataF t) {
    if (t < 0) {
        return 0.0;
    } else {
        return sin(2 * M_PI * omega * t);
    }
}

MyDataF CosineSource(MyDataF t) {
    if (t < 0) {
        return 0.0;
    } else {
        return cos(2 * M_PI * omega * t);
    }
}

MyDataF Source(MyDataF t) {
    if (pSource == NULL)return 0.0;
    else {
        return (*pSource)(t);
    }
}

MyDataF delays(MyDataF px, MyDataF py) {
    return ((px - xs) * cos_phi + (py - ys) * sin_phi);
}
