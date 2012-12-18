#include"common_data.h"


#define MDOUBLE MyDataF 
MDOUBLE *srcdat;
//c:\Documents and Settings\Wenbin Lin\My Documents\Visual Studio 2005\Projects\fdtd2d\fdtd2d\define_problem_space_and_parameters.h

void calsampos(int *xpos,int *ypos){
	MDOUBLE sample_x = 0;//where the sample field in x direction
	MDOUBLE sample_y = 0;//where the sample field in y direction
	
	sample_x = sample_y = lamda/10;

	if(IsTEx){
		*xpos = (int)(0.5+sample_x/dx)+nx/2;
		*ypos = (int)(0.5+sample_y/dy)+ny/2;
	}

	if(IsTMx){
		*xpos = (int)(0.5+sample_x/dx)+nx/2;
		*ypos = (int)(0.5+sample_y/dy)+ny/2;
	}
}

void InitCapture(int ttstep){
	int i;
	MyDataF t;
	//MyDataF *CapEf=NULL;
	FILE *src;
	MyDataF t_0;
	MyDataF tau = 1.5166/M_PI/f/3/5;
	t_0 =  1.8*tau;//dt*ttstep/10;
	omega=2*M_PI*f;
	
	printf("ttstep = %d\n",ttstep);
	if((src=fopen("source.dat","w"))==NULL)
	{
		fprintf(stderr,"Cannot open file source.dat!\n");
		exit(EXIT_FAILURE);
	}
	if((srcdat=(MyDataF*)malloc(ttstep*sizeof(MyDataF)))==NULL)
	{
		fprintf(stderr,"Cannot create space for srcdat!\n");
		exit(EXIT_FAILURE);
	}
	
	for(i=0,t=0;i<ttstep;i++,t+=dt){
		srcdat[i]=cos(omega*t);//exp(-pow((t-t_0)/tau,2));//
		//if(i<7800)
		//else srcdat[i]=0;
		fprintf(src,"%2.3e\t",srcdat[i]);
	}
	fclose(src);
}

void CaptureE(MyDataF *CapEf,int CurTimeStep,int i, int j){
	if(IsTEx){
		CapEf[CurTimeStep]=Ez_s.data[i*Ez_s.ny+j];
	}else if(IsTMx){
		CapEf[CurTimeStep]=Ey_s.data[i*Ey_s.ny+j];
	}
}
void scpfld(MyDataF *CapEf,int ttstep){
	FILE *fp=NULL;
	int i;

	if((fp=fopen("CapEfield.dat","w"))==NULL){
		printf("Cannot open field CapEfield.dat!\n");
		exit(0);
	}
	for(i=0;i<ttstep;i++)
		fprintf(fp,"%f\t",CapEf[i]);
	fprintf(fp,"\n");
	fclose(fp);

}
