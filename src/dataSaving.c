
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "commonData.h"
#include "fdtd.h"

MyDataF *srcdat = NULL;

void calsampos(int *xpos, int *ypos) {
    MyDataF sample_x = 0; //where the sample field in x direction
    MyDataF sample_y = 0; //where the sample field in y direction

    sample_x = sample_y = lamda / 10;

    if (IsTEx) {
        *xpos = (int) (0.5 + sample_x / dx) + nx / 2;
        *ypos = (int) (0.5 + sample_y / dy) + ny / 2;
    }

    if (IsTMx) {
        *xpos = (int) (0.5 + sample_x / dx) + nx / 2;
        *ypos = (int) (0.5 + sample_y / dy) + ny / 2;
    }
}

void InitCapture(int ttstep) {
    int i;
    MyDataF t;
    //MyDataF *CapEf=NULL;
    FILE *src = NULL;
    //MyDataF t_0;
    //MyDataF tau = 1.5166 / M_PI / f / 3 / 5;
    //t_0 = 1.8 * tau; //dt*ttstep/10;
    //omega = 2 * M_PI*f;

    printf("ttstep = %d\n", ttstep);
    if ((src = fopen("source.dat", "w")) == NULL) {
        fprintf(stderr, "Cannot open file source.dat!\n");
        exit(EXIT_FAILURE);
    }
    if ((srcdat = (MyDataF*) malloc(ttstep * sizeof (MyDataF))) == NULL) {
        fprintf(stderr, "Cannot create space for srcdat!\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0, t = 0; i < ttstep; i++, t += dt) {
        srcdat[i] = cos(omega * t); //exp(-pow((t-t_0)/tau,2));//
        //if(i<7800)
        //else srcdat[i]=0;
        fprintf(src, "%2.3e\t", srcdat[i]);
    }
    fclose(src);
}

void CaptureE(MyDataF *CapEf, int CurTimeStep, int i, int j) {
    if (IsTEx) {
        CapEf[CurTimeStep] = Ez_s.data[i * Ez_s.ny + j];
    } else if (IsTMx) {
        CapEf[CurTimeStep] = Ey_s.data[i * Ey_s.ny + j];
    }
}

void scpfld(MyDataF *CapEf, int ttstep) {
    FILE *fp = NULL;
    int i;

    if ((fp = fopen("CapEfield.dat", "w")) == NULL) {
        printf("Cannot open field CapEfield.dat!\n");
        exit(0);
    }
    for (i = 0; i < ttstep; i++) {
        fprintf(fp, "%f\t", CapEf[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);

}

void CheckIfOverLimit(int i, int j, MyDataF data) {
    //static int count = 0;
    if (fabs(data) > MAX_LIMIT) {

        printf("%d\t%d\t%-2.3e\n", i, j, data);
        //if(count%30 == 0){
        //////system("pause");
        //system("cls");
        //}
        //count++;
    }
}

void saveEfieldCenterTMz(int timestep) {
    static int cnt = 0;
    int jcenter = Hz_s.ny / 2, in;
    MyDataF ex, ey;

    char file[21] = "ec000000.dat";
    FILE *fp = NULL;
    if (timestep % SaveTimeStep != 0)
        return;
    file[2] = '0' + (cnt / 100000) % 10;
    file[3] = '0' + (cnt / 10000) % 10;
    file[4] = '0' + (cnt / 1000) % 10;
    file[5] = '0' + (cnt / 100) % 10;
    file[6] = '0' + (cnt / 10) % 10;
    file[7] = '0' + cnt % 10;
    if ((fp = fopen(file, "w")) == NULL) {
        printf("Cannot create file: \t%s\n", file);
        exit(EXIT_FAILURE);
    }
    for (in = tpis; in <= tpie && in < Hz_s.nx; in++) {
        ex = (Ex_s.data[in * Ex_s.ny + jcenter] + Ex_s.data[in * Ex_s.ny + jcenter + 1]) / 2;
        ey = (Ey_s.data[in * Ey_s.ny + jcenter] + Ey_s.data[(in + 1) * Ey_s.ny + jcenter]) / 2;
        fprintf(fp, "%4.4e\n", sqrt(ex * ex + ey * ey));
    }
    fclose(fp);
}

void CheckIfNeOverLimit(int i, int j, MyDataF data) {
    //static int count = 0;
    if (fabs(data) > NE_MAX_LIMIT) {

        printf("%d\t%d\t%-2.3e\n", i, j, data);
        //if(count%30 == 0){
        //////system("pause");
        //system("cls");
        //}
        //count++;
    }
}

void chkgrad(MyDataF str, MyDataF strpre, int i, int j) {

    MyDataF grad;

    grad = fabs(str - strpre);
    if (grad > 1e3) {
        printf("pre \t:%f\nnow\t:%f\ni = %d\tj = %d\ngrad = %f\n", strpre, str, i, j, grad);
        ////system("pause");
    }
}

void CheckIfNonZeros(MyStruct stru) {
    int i, j;
    for (i = 0; i < stru.nx; i++) {
        for (j = 0; j < stru.ny; j++) {
            if (fabs(stru.data[i * stru.ny + j]) > MAX_LIMIT) {
                printf("%d\t%d\t%-2.3e\n", i, j, stru.data[i * stru.ny + j]);
            }
        }
    }
}

//Save Fields

void SaveCapField(const int timestep) {
    static int cnt = 1;
    static int scnt = 1;
    int skipstep = 1e-9 / dt_F;
    //int leap = TotalTimeStep/SaveTimeStep;
    //int leap = 1;
    char file[21] = "ez000000.dat";

    if (scnt % skipstep == 0 /*&& cnt <=SAVE_LEAP*/) {

        file[2] = '0' + (cnt / 100000) % 10;
        file[3] = '0' + (cnt / 10000) % 10;
        file[4] = '0' + (cnt / 1000) % 10;
        file[5] = '0' + (cnt / 100) % 10;
        file[6] = '0' + (cnt / 10) % 10;
        file[7] = '0' + cnt % 10;
        //puts(fname);

        if (IsTEx) {
            CaptDataNoPML(cnt % 10, file, Ez_s, tpis, tpie, tpjs, tpje);
            //
            file[0] = 'e';
            file[1] = 'c';
            CaptDataCenter(cnt % 10, file, Ez_s, Ez_s.ny / 2, tpis, tpie);
            //file[0] = 'h';	file[1] = 'x';	CaptData(cnt%10,file,Hx_s); 
            //file[0] = 'h';	file[1] = 'y';	CaptData(cnt%10,file,Hy_s); 
        } else if (IsTMx) {
            //file[1] = 'x';	CaptData(cnt%10,file,Ex_s); 
            //file[1] = 'y';	CaptData(cnt%10,file,Ey_s);
            file[0] = 'h';
            file[1] = 'z';
            CaptDataNoPML(cnt % 10, file, Hz_s, tpis, tpie, tpjs, tpje);
            saveEfieldCenterTMz(timestep);
        }
        if(IfWithDensity) {
            file[0] = 'n';
            file[1] = 'e';
            CaptDataMNoPML(cnt % 10, file, ne, 2, tpis*m, tpie*m, tpjs*m, tpje * m);
            file[0] = 'n';
            file[1] = 'c';
            CaptDataCenter(cnt % 10, file, ne, ne.ny / 2, tpis*m, tpie * m);
            file[0] = 'e';
            file[1] = 'm';
            CaptDataMNoPML(cnt % 10, file, Erms, 2, tpis*m, tpie*m, tpjs*m, tpje * m);
            file[0] = 'm';
            file[1] = 'c';
            CaptDataCenter(cnt % 10, file, Erms, Erms.ny / 2, tpis*m, tpie * m);
        }
        cnt++;
    }
    scnt++;
}

void SaveErms(int nestep) {
    static int cnt = 1;
    static char ferms[21] = "er000000.dat";

    if (nestep % SAVE_ERMS_LEAP == 0) {
        ferms[2] = '0' + (cnt / 100000) % 10;
        ferms[3] = '0' + (cnt / 10000) % 10;
        ferms[4] = '0' + (cnt / 1000) % 10;
        ferms[5] = '0' + (cnt / 100) % 10;
        ferms[6] = '0' + (cnt / 10) % 10;
        ferms[7] = '0' + cnt % 10;
        //if(cnt!=30)CaptDataM(cnt%10,ferms,Erms,10);
        //else CaptData(cnt%10,ferms,Erms);
        CaptDataM(cnt % 10, ferms, Erms, 10);
        cnt++;
    }

}

void SaveRows(const MyStruct stru, int row, char *fname) {
    int i;
    FILE *fp;
    if (stru.data == NULL) {
        fprintf(stderr, "Null space to read.\n");
        printf("Save fails!\n");
        exit(EXIT_FAILURE);
    }
    if (row < 0 || row > stru.ny) {
        fprintf(stderr, "Invalid row.\n");
        printf("Save fails!\n");
        exit(EXIT_FAILURE);
    }
    if ((fp = fopen(fname, "w")) == NULL) {
        fprintf(stderr, "Cannot open file %s.\n", fname);
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < stru.ny; i++) {
        fprintf(fp, "%6.5e\t", stru.data[i * stru.ny + row]);
    }
    fprintf(fp, "\n");
    printf("Save file %s success.\n", fname);
    fclose(fp);
}


#define GetD(k,i,j) k.data[(i)*k.ny+j]


/******************************************************/

void SaveData(MyStruct data, char* filename, FILE *hfile);

/******************************************************/
void SaveToFile() {

    FILE *hfile;
    char dlist[15] = "DataList.dat";

    if ((hfile = fopen(dlist, "w")) == NULL) {
        printf("Cannot open file: \t%s\n", dlist);
        exit(0);
    }
    fprintf(hfile, "DataName\tnx\tny\n");

    if (IsTMx) {
        SaveData(Ex_s, "Ex_s.dat", hfile);
        SaveData(Ey_s, "Ey_s.dat", hfile);
        SaveData(Hz_s, "Hz_s.dat", hfile);

        SaveData(Vex, "Vex.dat", hfile);
        SaveData(Vey, "Vey.dat", hfile);
    }
    if (IsTEx) {
        SaveData(Hx_s, "Hx_s.dat", hfile);
        SaveData(Hy_s, "Hy_s.dat", hfile);
        SaveData(Ez_s, "Ez_s.dat", hfile);

        SaveData(Vez, "Vez.dat", hfile);
    }
    SaveData(ne, "ne.dat", hfile);
    fclose(hfile);
}

void SaveData(MyStruct data, char* filename, FILE *hfile) {
    int i, j;
    FILE *fp;

    if ((fp = fopen(filename, "w")) == NULL) {
        printf("Cannot create file: \t%s\n", filename);
        exit(0);
    }

    fprintf(hfile, "%s\t%d\t%d\n", filename, data.nx, data.ny);

    printf("\n%s\n", filename);
    //PrintData(data);
    //fwrite(data.data,sizeof(MyDataF),data.nx*data.ny,fp);
    for (j = 0; j < data.ny; j++) {
        for (i = 0; i < data.nx; i++)
            fprintf(fp, "%2.3E ", GetD(data, i, j));
        fprintf(fp, "\n");
    }

    fclose(fp);
}

#undef GetD

