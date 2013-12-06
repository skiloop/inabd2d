

#define SD(l,i,j) l.data[(i)*l.ny+j]
#define DD(s,p,i,j) (s.data[(i)*s.ny+j]+p.data[(i)*p.ny+j])


int Step = 1;
int MTimeStep = 1, MultiSize;

#include"updatescafield.h"
#include"captureefield.h"
#include"updatedensity.h"
#include"updatevelocity.h"
#include"connecting.h"
#include "updbnd.h"
//#include"updateincfield.h"

////////////////////////////////////////////////
///FUNCTION PREDEFINE
///////////////////////////////////////////////


void AdBound();
void InitSources();


/////////////////////////////////////////////
///END OF FUNCTION PREDEFING
/////////////////////////////////////////////

void scatfdtd() {

    int CurTimeStep; //current time step
    int xpos, ypos, sxpos, sypos;
    //FILE *fne;
    //int i,j;
    //MyDataF *CapEF;
    void PrintSourceSize(int, int, int, int);
    MyDataF cdtF;
    //int MTimeStep=1,MultiSize;
    //int cnt=0;
    int step_per_half_ns8 = (int) (0.5 + 0.125e-9 / dt_F);
    clock_t t_start, t_end;
    //MyDataF RealEz,rhtb;
    sxpos = (int) (0.5 + ny / 2); //x position of sources
    sypos = (int) (0.5 + ny / 2);
    //fne=fopen("cnep.dat","w");
    //y position of sources
    cnepx = m * ((int) (0.5 + 2.25 * lamda / dx) + tpis);
    cnepy = m * (nyp1 / 2);
    xpos = (int) (0.5 + nx / 2 + lamda / dx); //x position of field to be captured
    ypos = (int) (0.5 + ny / 2 + lamda / dy); //y position of field to be captured

    PrintSourceSize(xpos, ypos, sxpos, sypos);
    MultiSize = (int) (0.5 + dt_F / dt);

    printf("\ncnepx = %d\tcnepy = %d\n\n", cnepx, cnepy);
    printf("\nm\t=\t%d\nMultiSize\t=\t%d", m, MultiSize);
    printf("\nNumber of time steps: %d\n", TotalTimeStep);
    printf("Density Time Step: %d\n", Density_Time_Step * m);
    printf("Step per half ns: %d\n", step_per_half_ns8);
    printf("\n*************************************************************\n");
    //system("pause");

    //InitCapture(((TotalTimeStep>TotalTimeStep*MultiSize)?(TotalTimeStep*MultiSize):TotalTimeStep));
    //calsampos(&xpos,&ypos);
    //CapEF=(MyDataF*)malloc(TotalTimeStep*MultiSize*sizeof(MyDataF));
    //if(CapEF==NULL){
    //	printf("Cannot callocate space for CapEf!\n");
    //	exit(EXIT_FAILURE);
    //}

    sxpos = (int) (0.5 + nx / 2 + 0.125 * lamda / dx); //tpis+3;//tpis-SCATTER_FIELD_DOMAIN_BND_SIZE/2;//x position of sources
    sypos = (int) (0.5 + (tpjs + tpje) / 2); //(ny/2); //y position of sources
    filedat = fopen("neiter.dat", "w");
    if (filedat == NULL) {
        fprintf(stderr, "Cannot open neiter.dat for output!\n");
        exit(EXIT_FAILURE);
    }
    CurTime = -half_dt;
    InitSim();
    //InitIncFdtd();
    if (isConnect)initconnect(); ////InitFields();////Init fields before marching  loop i.e. at t = 0
    t_start = clock();
    for (CurTimeStep = 1; CurTimeStep <= Density_Time_Step && TotalTimeStep >= Step; CurTimeStep++) {

        for (cdtF = 0, MTimeStep = 1; MTimeStep <= MultiSize && cdtF < dt_F && TotalTimeStep >= Step; cdtF = cdtF + dt, MTimeStep++, Step++) {
            //MyDataF ezl,ezlr;

            CurTime += half_dt;

            UpdateMField();
            //connecting interface
            if (isConnect)mconnect(CurTime);
            //DispEMFields(CurTimeStep);
            //Hz_s.data[Hz_s.ny*Hz_s.ny/2+Hz_s.ny/2] += Source(CurTime);
            //UpdMFieldForPML();
            UpdateMFieldForPML(Hx_s, Hy_s, Hz_s, Ex_s, Ey_s, Ez_s);
            //DispEMFields(CurTimeStep);

            CurTime += half_dt;
            //DispEMFields(CurTimeStep);
            UpdateEField();
            if (!isConnect) {
                if (IsTEx)Ez_s.data[100 * Ez_s.ny + 50] += Source(CurTime);
                if (IsTMx)Hz_s.data[(pis + 10) * Hz_s.ny + 100] = Source(CurTime);
            }
            //DispEMFields(CurTimeStep);
            //connecting interface
            if (isConnect)econnect(CurTime);
            //Adjust_E_Field(CurTime-half_dt);////AddEInc(CurTime);//
            //DispEMFields(CurTimeStep);
            //UpdEFieldForPML();
            UpdateEFieldForPML(Hx_s, Hy_s, Hz_s, Ex_s, Ey_s, Ez_s);
            //DispEMFields(CurTimeStep);
            if (IfWithDensity) {
                UpdateVelocity();
                //CaptureE(CapEF,Step,xpos,ypos);
                if (if_erms_E_max == 1) {
                    CalSumESqrt_Emax();
                } else {
                    CalSumESqrt();
                }
            }
            if (Step == 100) {
                DispEMFields(CurTimeStep);
            }
        }
        if (IfWithDensity) {
            if (if_erms_E_max != 1) {
                CalErmsAtCoarseGrid();
            }
            InterpolatErms();
            SaveCapField(Step); // 保存数据
            UpdateDensity();
            UpdateCoeff();
            //SaveErms(CurTimeStep);
            //printf("ne = %15.14e\n",ne.data[cnepx*ne.ny+cnepy]);
            DispNe(CurTimeStep);
            //DispEMFields(CurTimeStep);
            ResetStructData(Erms);
        }
        printf("%d\t%f ns\n", Step, CurTime / 1e-9);
        ////
        //if(CurTimeStep%1==0)
    }//END FOR
    //SaveD(cnt);
    t_end = clock();
    printf("Total time used : %ld\n", t_end - t_start);
    //system("pause");
    EndMatlabSim();
    ////Save Sample Field
    //scpfld(CapEF,TotalTimeStep);
    //free(CapEF);
    end_connect();
    free(srcdat);
    fclose(filedat);

}//END OF FDTDLOOP

void PrintSourceSize(int xpos, int ypos, int sxpos, int sypos) {
    printf("\nnx\t=\t%d\nny\t=\t%d\nxpos\t=\t%d\nypos\t=\t%d\nsxpos\t=\t%d\nsypos\t=\t%d\nnbound\t=\t%d\n", nx, ny, xpos, ypos, sxpos, sypos, nbound);
}
#undef SD
