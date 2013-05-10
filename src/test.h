#include "common_data.h"

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
    for (i = 0; i < stru.nx; i++)
        for (j = 0; j < stru.ny; j++)
            if (fabs(stru.data[i * stru.ny + j]) > MAX_LIMIT)
                printf("%d\t%d\t%-2.3e\n", i, j, stru.data[i * stru.ny + j]);
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
        file[0] = 'n';
        file[1] = 'e';
        CaptDataMNoPML(cnt % 10, file, ne, 4, tpis*m, tpie*m, tpjs*m, tpje * m);
        file[0] = 'n';
        file[1] = 'c';
        CaptDataCenter(cnt % 10, file, ne, ne.ny / 2, tpis*m, tpie * m);
        file[0] = 'e';
        file[1] = 'm';
        CaptDataMNoPML(cnt % 10, file, Erms, 4, tpis*m, tpie*m, tpjs*m, tpje * m);
        file[0] = 'm';
        file[1] = 'c';
        CaptDataCenter(cnt % 10, file, Erms, Erms.ny / 2, tpis*m, tpie * m);

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
    for (i = 0; i < stru.ny; i++)
        fprintf(fp, "%6.5e\t", stru.data[i * stru.ny + row]);
    fprintf(fp, "\n");
    printf("Save file %s success.\n", fname);
    fclose(fp);
}


//MyDataF vez[TotalTimeStep],ez[TotalTimeStep],cne[TotalTimeStep];
//MyDataF vez2[TotalTimeStep],ez2[TotalTimeStep],cne2[TotalTimeStep];
//void Cap(int cxpos,int cypos,int cnt){
//
//	vez[cnt] =Vez.data[cxpos*Vez.ny+cypos];
//	ez[cnt] = Ez.data[cxpos*Vez.ny+cypos];
//	cne[cnt] = ne.data[m*cxpos*ne.ny+m*cypos];
//	vez2[cnt] =Vez.data[(cxpos)*Vez.ny+cypos+45];
//	ez2[cnt] = Ez.data[cxpos*Vez.ny+cypos+45];
//	cne2[cnt] = ne.data[m*cxpos*ne.ny+m*cypos+45*m];
//}
//void SaveD(int cnt){
//	FILE *fp;
//	int i;
//	fp=fopen("cvez.dat","w");
//	for(i=0;i<cnt;i++)
//		fprintf(fp,"%5.2e\t",vez[i]);
//	fclose(fp);
//
//	fp=fopen("cez.dat","w");
//	for(i=0;i<cnt;i++)
//		fprintf(fp,"%5.2e\t",ez[i]);
//	fclose(fp);
//	fp=fopen("cne.dat","w");
//	for(i=0;i<cnt;i++)
//		fprintf(fp,"%5.2e\t",cne[i]);
//	fclose(fp);
//
//	fp=fopen("cvez2.dat","w");
//	for(i=0;i<cnt;i++)
//		fprintf(fp,"%5.2e\t",vez2[i]);
//	fclose(fp);
//
//	fp=fopen("cez2.dat","w");
//	for(i=0;i<cnt;i++)
//		fprintf(fp,"%5.2e\t",ez2[i]);
//	fclose(fp);
//	fp=fopen("cne2.dat","w");
//	for(i=0;i<cnt;i++)
//		fprintf(fp,"%5.2e\t",cne2[i]);
//	fclose(fp);
//}
