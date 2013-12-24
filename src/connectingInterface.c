
#include <stdlib.h>
#include <math.h>

#include "common.h"
#include "fdtd.h"
#include "dataType.h"
#include "commonData.h"
#include "connectingInterface.h"


MyDataF *Ei, *Hi;
int eilen, hilen;
MyDataF CMur;
int start_index_x;
int start_index_y;

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

MyDataF cos_phi, sin_phi;
MyDataF t_per_cell; //time wave propagates per Yee cell
int xs, ys, xe, ye; //

#define INC_SIZE 0

void initSource() {
    switch (srcType) {
        case GAUSSIAN:
            tauSquare = T * T / 4;
            t_0 = T *1.2;
            pSource = GaussianSource;
            break;
        case COSINE:
            pSource=CosineSource;
            break;
        case SINE:
        default:
            pSource=SineSource;
    }
}

void initconnect() {
    //#ifdef  MATLAB_SIMULATION
    //	mxArray *MyArray;
    //
    //#endif
    int i;
    MyDataF ds;
    MyDataF PhaseVelRatio(MyDataF);

    ds = dx;
    eilen = tpie - tpis + 5;
    hilen = eilen - 1;
    start_index_x = -tpis;
    start_index_y = -tpjs;

    CMur = (c * dt - ds) / (c * dt + ds);
    if (IsTMx) {
        Chzey = (Ratio_y) * dt / mu_0 / dx;
        Chzex = (Ratio_x) * dt / mu_0 / dy;
        Ceyhz = -dt / ds / eps_0;
        Cexhz = dt / ds / eps_0;
        chiei = -dt / mu_0 / dx;
        ceihi = -dt / eps_0 / dx;
    }
    if (IsTEx) {
        chiei = dt / mu_0 / dx;
        ceihi = dt / eps_0 / dx;
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
    //	MyArray=mxCreateDoubleMatrix(1,eilen,mxREAL);
    //	memcpy(mxGetPr(MyArray), Ei,eilen*sizeof(MyDataF));
    //	engPutVariable(ep,"Ei",MyArray);
    //	engEvalString(ep,"figure,eih=plot(Ei);grid on;set(gca,'ylim',[-10 10]);");//%surf(data1);");//
    //	
    //	mxDestroyArray(MyArray);
    //#endif

    //calculate phase velocity ratio
    Vp_ratio = PhaseVelRatio(0) / PhaseVelRatio(phi);

}

void econnect(MyDataF t) {
    //#ifdef  MATLAB_SIMULATION
    //	mxArray *MyArray;
    //
    //#endif
    int i, index, ind;
    int di;

    MyDataF ei_last, ei_last2, df, d;
    //MyDataF ds;

    //ds = dx;
    ei_last = Ei[eilen - 1];
    ei_last2 = Ei[eilen - 2];
    Ei[0] = E0 * Source(t);
    for (i = 1; i <= eilen - 2; i++)
        Ei[i] = Ei[i] + ceihi * (Hi[i] - Hi[i - 1]);

    //mur boundary 
    Ei[eilen - 1] = ei_last2 + CMur * (Ei[eilen - 2] - ei_last);


    //#ifdef  MATLAB_SIMULATION
    //	MyArray=mxCreateDoubleMatrix(1,eilen,mxREAL);
    //	memcpy(mxGetPr(MyArray), Ei,eilen*sizeof(MyDataF));
    //	engPutVariable(ep,"Ei",MyArray);
    //	engEvalString(ep,"set(eih,'ydata',Ei);");//surf(data1);");//
    //	//engEvalString(ep,"title(ind);ind=ind+1;");
    //	mxDestroyArray(MyArray);
    //#endif
    if (IsTEx) {
        for (ind = tpis; ind <= tpie; ind++) {
            index = ind * Ez.ny;
            //bottom 
            d = Dhxb[ind + start_index_x] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ez.data[index + tpjs] += Cehz.data[index + tpjs] * Ratio_x * (Hi[di] + df * (Hi[di + 1] - Hi[di]));

            //top
            d = Dhxt[ind + start_index_x] + 0.5;

            di = (int) floor(d);
            df = d - di;
            Ez.data[index + tpje] -= Cehz.data[index + tpje] * Ratio_x * (Hi[di] + df * (Hi[di + 1] - Hi[di]));
        }
        //bound xn,xp
        for (ind = tpjs; ind <= tpje; ind++) {
            //left side
            index = tpis * Ez.ny + ind;
            d = Dhyl[ind + start_index_y] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ez.data[index] -= Cehz.data[index]*(-Ratio_y)*(Hi[di] + df * (Hi[di + 1] - Hi[di]));

            //right side
            index = tpie * Ez.ny + ind;
            d = Dhyr[ind + start_index_y] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ez.data[index] += Cehz.data[index]*(-Ratio_y)*(Hi[di] + df * (Hi[di + 1] - Hi[di]));
        }
    }
    if (IsTMx) {
        int i, j;
        int ind;
        for (ind = 0, j = tpjs; j < tpje; j++, ind++) {
            d = Dhzl[ind] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ey.data[tpis * Ey.ny + j] -= Ceyhz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //left side

            d = Dhzr[ind] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ey.data[tpie * Ey.ny + j] += Ceyhz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //right side

        }
        for (ind = 0, i = tpis; i < tpie; i++, ind++) {
            d = Dhzb[ind] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ex.data[i * Ex.ny + tpjs] -= Cexhz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //bottom side

            d = Dhzt[ind] + 0.5;
            di = (int) floor(d);
            df = d - di;
            Ex.data[i * Ex.ny + tpje] += Cexhz * (Hi[di] + df * (Hi[di + 1] - Hi[di])); //top side
        }
    }
}

void mconnect(MyDataF t) {

    int i, ind, ind1, ind3;
    MyDataF hi_last, hi_last2;
    //MyDataF ds;

    int di;
    MyDataF df;
    //MyDataF tmp1;

    //ds = dx;
    hi_last = Hi[hilen - 1];
    hi_last2 = Hi[hilen - 2];
    for (i = 0; i < hilen - 1; i++)
        Hi[i] = Hi[i] + chiei * (Ei[i + 1] - Ei[i]);


    //Mur boundary
    Hi[hilen - 1] = hi_last2 + CMur * (Hi[hilen - 2] - hi_last);
    if (IsTMx) {
        int mtpis, mtpjs;
        mtpis = tpis - 1;
        mtpjs = tpjs - 1;
        for (ind = tpis, ind1 = 0, ind3 = tpje; ind < tpie; ind++, ind1++) {
            di = (unsigned) floor(Dexb[ind1]);
            df = Dexb[ind1] - di;
            Hz.data[ind * Hz.ny + mtpjs] -= Chzex * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //bottom side
            di = (unsigned) floor(Dext[ind1]);
            df = Dext[ind1] - di;
            Hz.data[ind * Hz.ny + ind3] += Chzex * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //top side
        }
        for (ind = tpjs, ind1 = 0, ind3 = tpie; ind < tpje; ind++, ind1++) {
            di = (unsigned) floor(Deyl[ind1]);
            df = Deyl[ind1] - di;
            Hz.data[mtpis * Hz.ny + ind] -= Chzey * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //left side
            di = (unsigned) floor(Deyr[ind1]);
            df = Deyr[ind1] - di;
            Hz.data[ind3 * Hz.ny + ind] += Chzey * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //right side
        }

    }
    if (IsTEx) {
        for (ind = tpjs; ind <= tpje; ind++) {
            //left side
            di = (int) floor(Dezl[ind + start_index_y]);
            df = Dezl[ind + start_index_y] - di;
            di = di + 1;
            Hy.data[(tpis - 1) * Hy.ny + ind] -= Chyez * (Ei[di] + df * (Ei[di + 1] - Ei[di])); //0;//
            //right side
            di = (int) floor(Dezr[ind + start_index_y]);
            df = Dezr[ind + start_index_y] - di;
            di = di + 1;
            Hy.data[tpie * Hy.ny + ind] += Chyez * (Ei[di] + df * (Ei[di + 1] - Ei[di]));
        }
        //Adjust Hy
        for (ind = tpis; ind <= tpie; ind++) {
            // bottom
            di = (int) floor(Dezb[ind + start_index_x]);
            df = Dezb[ind + start_index_x] - di;
            di = di + 1;

            Hx.data[ind * Hx.ny + tpjs - 1] -= Chxez * (Ei[di] + df * (Ei[di + 1] - Ei[di]));
            //yp
            di = (int) floor(Dezt[ind + start_index_x]);
            df = Dezt[ind + start_index_x] - di;
            di = di + 1;
            Hx.data[ind * Hx.ny + tpje] += Chxez * (Ei[di] + df * (Ei[di + 1] - Ei[di]));
        }
    }

}

void end_connect() {
    free(Ei);
    free(Hi);
}

MyDataF PhaseVelRatio(MyDataF angle) {
    MyDataF A, B, C, S, k, N, kp;

    N = lamda / dx;
    S = c * dt / dx;
    A = 0.5 * dx * cos(angle);
    B = 0.5 * dx * sin(angle);
    C = sin(M_PI * S / N) * sin(M_PI * S / N) / S / S;

    k = 0.5;
    kp = 1;
    while (fabs(k - kp) > 1e-4) {
        kp = k;
        k = kp - (sin(A * kp) * sin(A * kp) + sin(B * kp) * sin(B * kp) - C) / (A * sin(2 * A * kp) + B * sin(2 * B * kp));
    }

    return 2 * M_PI * c / k;
}

void CalDelays() {
    int i;
    int Max_Index_Dez_x;
    int Max_Index_Dez_y;

    int Max_Index_Dhx;
    int Max_Index_Dhy;

    //FILE *fp1,*fp2;
    MyDataF delays(MyDataF px, MyDataF py);


    t_per_cell = dx / c;
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    phi = phi - 2 * M_PI * (floor(phi / (2 * M_PI)));
    if ((phi <= 0.5 * M_PI)) {
        xs = tpis;
        ys = tpjs;
        xe = tpie;
        ye = tpje;
    } else if ((phi <= M_PI)) {
        xs = tpie;
        ys = tpjs;
        xe = tpis;
        ye = tpje;
    } else if ((phi <= 1.5 * M_PI)) {
        xs = tpie;
        ys = tpje;
        xe = tpis;
        ye = tpjs;
    } else {
        xs = tpis;
        ys = tpje;
        xe = tpie;
        ye = tpjs;
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
            Dhxb[i] = delays(tpis + i, tpjs - 0.5);
            Dhxt[i] = delays(tpis + i, tpje + 0.5);
        }
        for (i = 0; i <= Max_Index_Dhy; i++) {
            Dhyl[i] = delays(tpis - 0.5, tpjs + i);
            Dhyr[i] = delays(tpie + 0.5, tpjs + i);
        }
        for (i = 0; i <= Max_Index_Dez_x; i++) {
            Dezb[i] = delays(tpis + i, tpjs);
            Dezt[i] = delays(tpis + i, tpje);
        }
        for (i = 0; i <= Max_Index_Dez_y; i++) {
            Dezl[i] = delays(tpis, tpjs + i);
            Dezr[i] = delays(tpie, tpjs + i);
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
            Dexb[i] = delays(tpis + i + 0.5, tpjs) + 1; //plus 1 is because orignal point is the second one i.e. Ei[1]
            Dext[i] = delays(tpis + i + 0.5, tpje) + 1; //plus 1 is because orignal point is the second one i.e. Ei[1]
        }
        for (i = 0; i < Max_Index_Dey; i++) {
            Deyl[i] = delays(tpis, tpjs + i + 0.5) + 1; //plus 1 is because orignal point is the second one i.e. Ei[1]
            Deyr[i] = delays(tpie, tpjs + i + 0.5) + 1; //plus 1 is because orignal point is the second one i.e. Ei[1]
        }
        for (i = 0; i < Max_Index_Dhz_x; i++) {
            Dhzb[i] = delays(tpis + i + 0.5, tpjs - 0.5);
            Dhzt[i] = delays(tpis + i + 0.5, tpje + 0.5);
        }
        for (i = 0; i < Max_Index_Dhz_y; i++) {
            Dhzl[i] = delays(tpis - 0.5, tpjs + i + 0.5);
            Dhzr[i] = delays(tpie + 0.5, tpjs + i + 0.5);
        }
    }

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

MyDataF Source(MyDataF t){
    if (pSource==NULL)return 0.0;
    else {
        return (*pSource)(t);
    }
}

void InitFields() {
    //Init Fields at t = 0;
    int j;
    for (j = tpjs; j <= tpje; j++)
        Ez.data[tpis * Ez.ny + j] = Ez0 * sin(omega * CurTime);
}
//

MyDataF delays(MyDataF px, MyDataF py) {
    return ((px - xs) * cos_phi + (py - ys) * sin_phi);
}

void AddEInc(MyDataF t) {
    int i;
    for (i = tpis; i <= tpie; i++) {
        Ez.data[i * Ez.ny + tpjs] = E0 * Source(t - Dezb[i - tpis] * t_per_cell);
        Ez.data[i * Ez.ny + tpje] = E0 * Source(t - Dezt[i - tpis] * t_per_cell);
    }
    for (i = tpjs + 1; i < tpje; i++) {
        Ez.data[tpis * Ez.ny + i] = E0 * Source(t - Dezl[i - tpjs] * t_per_cell);
        Ez.data[tpie * Ez.ny + i] = E0 * Source(t - Dezr[i - tpjs] * t_per_cell);
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

void IntTtlFldDmnBnd() {
    //this function define where the total field domain is
    //Note: pis,pjs,pie,pje must be defined and initialzed 
    //		before calling this function
    tpis = pis + SCATTER_FIELD_DOMAIN_BND_SIZE;
    tpjs = pjs + SCATTER_FIELD_DOMAIN_BND_SIZE;
    tpie = pie - SCATTER_FIELD_DOMAIN_BND_SIZE;
    tpje = pje - SCATTER_FIELD_DOMAIN_BND_SIZE;
    printf("\n*************************************************************\n");
    printf("\ntpis\t=\t%d\ntpjs\t=\t%d\ntpie\t=\t%d\ntpje\t=\t%d\n", tpis, tpjs, tpie, tpje);
    CalDelays();
}
//
//void Adjust_E_Field(MyDataF pre_t) {//adjust TEx electric field which at total fields boundary
//    int ind;
//    int index;
//    //bound yn,yp
//
//    for (ind = tpis; ind <= tpie; ind++) {
//        index = ind * Ez_s.ny;
//        //bottom  Cehz.data[index+tpjs]
//        Ez_s.data[index + tpjs] += 0; //Cehz.data[index+tpjs]/dy*Hx0*Source(pre_t-Dhxb[ind-tpis]*t_per_cell);
//        //top
//        Ez_s.data[index + tpje] -= 0; //Cehz.data[index+tpje]/dy*Hx0*Source(pre_t-Dhxt[ind-tpis]*t_per_cell);
//    }
//    //bound xn,xp Cehz.data[index+tpje]
//    for (ind = tpjs; ind <= tpje; ind++) {
//        index = tpis * Ez_s.ny + ind;
//        //left sideCehz.data[index]
//        Ez_s.data[index] -= Cehz.data[index] / dx * Hy0 * Source(pre_t - Dhyl[ind - tpjs] * t_per_cell);
//        index = tpie * Ez_s.ny + ind;
//        //right
//        Ez_s.data[index] += Cehz.data[index] / dx * Hy0 * Source(pre_t - Dhyr[ind - tpjs] * t_per_cell);
//    }
//
//}
//
////
//
//void Adjust_M_Field(MyDataF pre_t) {
//    //adjust magnetic field of TMx in total field boundary
//    int ind;
//    // MyDataF tmp;
//    //Adjust Hy
//    for (ind = tpjs; ind <= tpje; ind++) {
//        //left
//        Hy_s.data[(tpis - 1) * Hy_s.ny + ind] -= Chyez * Ez0 * Source(pre_t - Dezl[ind - tpjs] * t_per_cell);
//        //right
//        Hy_s.data[tpie * Hy_s.ny + ind] += Chyez * Ez0 * Source(pre_t - Dezr[ind - tpjs] * t_per_cell);
//    }
//    //Adjust Hy
//    for (ind = tpis; ind <= tpie; ind++) {
//        MyDataF Einc;
//
//        //bottom
//        Einc = Ez0 * Source(pre_t - Dezb[ind - tpis] * t_per_cell);
//        tmp = Chxez*Einc;
//        Hx_s.data[ind * Hx_s.ny + tpjs - 1] -= 0; //tmp;//0;
//        //top
//        Hx_s.data[ind * Hx_s.ny + tpje] += 0; //Chxez*Ez0*Source(pre_t-Dezt[ind-tpis]*t_per_cell);
//    }
//}
