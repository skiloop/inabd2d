#ifndef __SOURCE_H__
#define __SOURCE_H__
//#include <iostream>
//using namespace std;
//Delay of time comparing with position tpis
//'D' for delay;
//'ez' for electric field Ez;
//'hx' 'hy' for magnetic field Hx and Hy
//'t' for top;'b' for bottom; 'l' for left;'r' for right
//Ez Hx has delay arrays  

//Hy Delay
MyDataF *Dhyl,*Dhyr;
//Hx Delay
MyDataF *Dhxb,*Dhxt;
//Ez Delay
MyDataF *Dezb,*Dezt,*Dezl,*Dezr;

//Ey Delay
MyDataF *Deyl,*Deyr;
//Ex Delay
MyDataF *Dexb,*Dext;
//Hz Delay
MyDataF *Dhzb,*Dhzt,*Dhzl,*Dhzr;

MyDataF cos_phi,sin_phi;
MyDataF t_per_cell; //time wave propagates per Yee cell
int xs,ys,xe,ye;//


void CalDelays()
{
	int i;
	int Max_Index_Dez_x;
	int Max_Index_Dez_y;

	int Max_Index_Dhx;
	int Max_Index_Dhy;

	//FILE *fp1,*fp2;
	MyDataF delays(MyDataF px,MyDataF py);
	
	
	t_per_cell = dx/c;
	cos_phi=cos(phi);
	sin_phi=sin(phi);

	phi = phi-2*M_PI*(floor(phi/(2*M_PI)));
	if((phi<=0.5*M_PI)){
		xs = tpis;		ys = tpjs;
		xe = tpie;		ye = tpje;
	}else if((phi<=M_PI)){
		xs = tpie;		ys = tpjs;
		xe = tpis;		ye = tpje;
	}else if((phi<=1.5*M_PI)){
		xs = tpie;		ys = tpje;
		xe = tpis;		ye = tpjs;
	}else {
		xs = tpis;		ys = tpje;
		xe = tpie;		ye = tpjs;
	}
	if(IsTEx){
		Max_Index_Dhx = Max_Index_Dez_x = abs(xe-xs);
		Max_Index_Dhy = Max_Index_Dez_y = abs(ye-ys);

		Dhxb = (MyDataF*)malloc((Max_Index_Dhx+1)*sizeof(MyDataF));
		Dhxt = (MyDataF*)malloc((Max_Index_Dhx+1)*sizeof(MyDataF));
		if(Dhxb==NULL||Dhxt==NULL){
			printf("Cannot allocate enough memory for Dhx!\n");
			exit(EXIT_FAILURE);
		}
		Dhyl = (MyDataF*)malloc((Max_Index_Dhy+1)*sizeof(MyDataF));
		Dhyr = (MyDataF*)malloc((Max_Index_Dhy+1)*sizeof(MyDataF));
		if(Dhyl==NULL||Dhyr==NULL){
			printf("Cannot allocate enough memory for Dhy!\n");
			exit(EXIT_FAILURE);
		}

		Dezb = (MyDataF*)malloc((Max_Index_Dez_x+1)*sizeof(MyDataF));
		Dezt = (MyDataF*)malloc((Max_Index_Dez_x+1)*sizeof(MyDataF));
		Dezl = (MyDataF*)malloc((Max_Index_Dez_y+1)*sizeof(MyDataF));
		Dezr = (MyDataF*)malloc((Max_Index_Dez_y+1)*sizeof(MyDataF));

		if(Dezt==NULL||Dezb==NULL||Dezl==NULL||Dezr==NULL){
			printf("Cannot allocate enough memory for Dez!\n");
			exit(EXIT_FAILURE);
		}

		for(i=0;i<=Max_Index_Dhx;i++)
		{
			Dhxb[i]=delays(tpis+i,tpjs-0.5);
			Dhxt[i]=delays(tpis+i,tpje+0.5);
		}
		for(i=0;i<=Max_Index_Dhy;i++)
		{
			Dhyl[i]=delays(tpis-0.5,tpjs+i);
			Dhyr[i]=delays(tpie+0.5,tpjs+i);
		}
		for(i=0;i<=Max_Index_Dez_x;i++){
			Dezb[i]=delays(tpis+i,tpjs);
			Dezt[i]=delays(tpis+i,tpje);
		}
		for(i=0;i<=Max_Index_Dez_y;i++){
			Dezl[i]=delays(tpis,tpjs+i);
			Dezr[i]=delays(tpie,tpjs+i);
		}
	}
	if(IsTMx){
		//Ex Ey are total fields and Hz is scattered field
		int Max_Index_Dhz_x;
		int Max_Index_Dhz_y;

		int Max_Index_Dex;
		int Max_Index_Dey;

		Max_Index_Dex = Max_Index_Dhz_x = (int)abs((double)(xe-xs));
		Max_Index_Dey = Max_Index_Dhz_y = (int)abs((double)(ye-ys));

		Dexb = (MyDataF*)calloc((Max_Index_Dex),sizeof(MyDataF));
		Dext = (MyDataF*)calloc((Max_Index_Dex),sizeof(MyDataF));

		if(Dexb==NULL || Dext==NULL){
			printf("Cannot allocate enough memory for Dex!");
			exit(EXIT_FAILURE);
		}
		Deyl = (MyDataF*)calloc((Max_Index_Dey),sizeof(MyDataF));
		Deyr = (MyDataF*)calloc((Max_Index_Dey),sizeof(MyDataF));
		if(Deyl==NULL || Deyr==NULL){
			printf("Cannot allocate enough memory for Dey!\n");;
			exit(EXIT_FAILURE);
		}

		Dhzb = (MyDataF*)calloc((Max_Index_Dhz_x),sizeof(MyDataF));
		Dhzt = (MyDataF*)calloc((Max_Index_Dhz_x),sizeof(MyDataF));
		Dhzl = (MyDataF*)calloc((Max_Index_Dhz_y),sizeof(MyDataF));
		Dhzr = (MyDataF*)calloc((Max_Index_Dhz_y),sizeof(MyDataF));
		if(Dhzt==NULL||Dhzb==NULL||Dhzl==NULL||Dhzr==NULL){
			fprintf(stderr,"Cannot allocate enough memory for Dhz!");
			exit(EXIT_FAILURE);
		}

		for(i=0;i<Max_Index_Dex;i++)
		{
			Dexb[i] = delays(tpis+i+0.5, tpjs)+1;//plus 1 is because orignal point is the second one i.e. Ei[1]
			Dext[i] = delays(tpis+i+0.5, tpje)+1;//plus 1 is because orignal point is the second one i.e. Ei[1]
		}
		for(i=0;i<Max_Index_Dey;i++)
		{
			Deyl[i]=delays(tpis,tpjs+i+0.5)+1;//plus 1 is because orignal point is the second one i.e. Ei[1]
			Deyr[i]=delays(tpie,tpjs+i+0.5)+1;//plus 1 is because orignal point is the second one i.e. Ei[1]
		}
		for(i=0;i<Max_Index_Dhz_x;i++){
			Dhzb[i]=delays(tpis+i+0.5,tpjs-0.5);
			Dhzt[i]=delays(tpis+i+0.5,tpje+0.5);
		}
		for(i=0;i<Max_Index_Dhz_y;i++){
			Dhzl[i]=delays(tpis-0.5,tpjs+i+0.5);
			Dhzr[i]=delays(tpie+0.5,tpjs+i+0.5);
		}
	}
	////save delays
	//fp1=fopen("dhxb.dat","w");
	//fp2=fopen("dhxt.dat","w");
	//for(i=0;i<=Max_Index_Dhx;i++)
	//{
	//	fprintf(fp1,"%12.10e\t",Dhxb[i]);
	//	fprintf(fp2,"%12.10e\t",Dhxt[i]);
	//}
	//fclose(fp1);
	//fclose(fp2);


	//fp1=fopen("dhyl.dat","w");
	//fp2=fopen("dhyr.dat","w");
	//for(i=0;i<=Max_Index_Dhy;i++)
	//{
	//	fprintf(fp1,"%12.10e\t",Dhyl[i]);
	//	fprintf(fp2,"%12.10e\t",Dhyr[i]);
	//}
	//fclose(fp1);fclose(fp2);
	//
	//fp1=fopen("dezb.dat","w");
	//fp2=fopen("dezt.dat","w");
	//for(i=0;i<=Max_Index_Dez_x;i++){
	//	fprintf(fp1,"%12.10e\t",Dezb[i]);
	//	fprintf(fp2,"%12.10e\t",Dezt[i]);
	//}
	//fclose(fp1);fclose(fp2);
	//
	//
	//fp1=fopen("dezl.dat","w");
	//fp2=fopen("dezr.dat","w");
	//for(i=0;i<=Max_Index_Dez_y;i++){
	//	fprintf(fp1,"%12.10e\t",Dezl[i]);
	//	fprintf(fp2,"%12.10e\t",Dezr[i]);
	//}
	//fclose(fp1);fclose(fp2);
	

}
MyDataF Source(MyDataF t)
{
	MyDataF t_0;
	MyDataF tau = 3.1E-9;
	t_0 =  0.8*tau;//dt*ttstep/10;//||t>1000*dt
	//
	if(t<0)return 0;
	return cos(omega*t);//1;//-cos(2*M_PI*f*t)*exp(-4*M_PI*pow((t-t_0)/tau,2));//
	//return 0;1.0;
}

void InitFields()
{
	//Init Fields at t = 0;
	int j;
	for(j=tpjs;j<=tpje;j++)
		Ez_s.data[tpis*Ez_s.ny+j] = Ez0*sin(omega*CurTime);
}
//
MyDataF delays(MyDataF px,MyDataF py)
{
	return ((px-xs)*cos_phi+(py-ys)*sin_phi);
}
void AddEInc(MyDataF t)
{
	int i;
	for(i=tpis;i<=tpie;i++)
	{
		Ez_s.data[i*Ez_s.ny+tpjs] = E0*Source(t-Dezb[i-tpis]*t_per_cell);
		Ez_s.data[i*Ez_s.ny+tpje] = E0*Source(t-Dezt[i-tpis]*t_per_cell);
	}
	for(i=tpjs+1;i<tpje;i++)
	{
		Ez_s.data[tpis*Ez_s.ny+i] = E0*Source(t-Dezl[i-tpjs]*t_per_cell);
		Ez_s.data[tpie*Ez_s.ny+i] = E0*Source(t-Dezr[i-tpjs]*t_per_cell);
	}
}
void FreeDelayArrays()
{
	if(IsTEx){
		free(Dezb);
		free(Dezt);
		free(Dezl);
		free(Dezr);

		free(Dhxb);
		free(Dhxt);
		free(Dhyl);
		free(Dhyr);
	}
	if(IsTMx)
	{
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

#endif 
