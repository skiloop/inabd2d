#ifndef CALCULATE_DOMAIN_SIZE
#define CALCULATE_DOMAIN_SIZE

#include"common_data.h"
//#include"init_com_data.h"

void CalDomainSize(){

	int nxi,nyi,nxip1,nyip1;;
	nxp1=nx=nxm1=ny=nyp1=nym1=0;

	nxi = (int)(0.5+(int)(0.5+lamda /dx) * NUMBER_OF_WAVELENGTHS_IN_DOMAIN);
	nyi = (int)(0.5+(int)(0.5+lamda /dy) * NUMBER_OF_WAVELENGTHS_IN_DOMAIN);
	nxip1 = nxi + 1;
	nyip1 = nyi + 1;

	nx= nxi  + 2*(NUMBER_OF_CELLS_IN_PML_BOUND + SCATTER_FIELD_DOMAIN_BND_SIZE);
	ny= nyi  + 2*(NUMBER_OF_CELLS_IN_PML_BOUND + SCATTER_FIELD_DOMAIN_BND_SIZE);


	nxp1=nx+1;
	nyp1=ny+1;
	nxm1=nx-1;
	nym1=ny-1;

	

	if(IsTMx){

		InitMyStr(nx,nyp1,&Ex_s_pre);	
		//InitMyStr(nx,ny,&Hz_s_pre);
		InitMyStr(nxp1,ny,&Ey_s_pre);

		//InitMyStr(nx,nyp1,&Ex);
		//InitMyStr(nxp1,ny,&Ey);
		//InitMyStr(nx,ny,&Hz);

		InitMyStr(nx,nyp1,&Ex_i);
		InitMyStr(nxp1,ny,&Ey_i);
		InitMyStr(nx,ny,&Hz_i);
		/*----------------------------------*/
		
		InitMyStr(nx,nyp1,&Ex_s);
	
		//InitMyStr(nx,nyp1,&Ex_i);
		//InitMyStr(nx,nyp1,&Ex_i_pre);
		InitMyStr(nx,nyp1,&Vex);
		InitMyStr(nx,nyp1,&Cevx);
		InitMyStr(nx,nyp1,&Ceex);
		InitMyStr(nx,nyp1,&Cehx);
		/*----------------------------------*/


		/*----------------------------------*/

		
		InitMyStr(nxp1,ny,&Ey_s);
		//InitMyStr(nxp1,ny,&Ey_i);
		//InitMyStr(nxp1,ny,&Ey_i_pre);	

		InitMyStr(nxp1,ny,&Vey);
		InitMyStr(nxp1,ny,&Cevy);
		InitMyStr(nxp1,ny,&Ceey);
		InitMyStr(nxp1,ny,&Cehy);
		/*----------------------------------*/


		/*----------------------------------*/
		
		InitMyStr(nx,ny,&Hz_s);
		//InitMyStr(nx,ny,&Hz_i);
		//InitMyStr(nx,ny,&Hz_i_pre);
		/*----------------------------------*/
	}
	if(IsTEx){

		
		//InitMyStr(nxi,nyi,&Ez_i);
		//InitMyStr(nxi,nyip1,&Hx_i);
		//InitMyStr(nxip1,nyi,&Hx_i);

		//InitMyStr(nxp1,ny,&Hx_s_pre);
		//InitMyStr(nx,nyp1,&Hy_s_pre);
		InitMyStr(nxp1,nyp1,&Ez_s_pre);



		InitMyStr(nxp1,ny,&Hx_i);
		InitMyStr(nx,nyp1,&Hy_i);
		InitMyStr(nxp1,nyp1,&Ez_i);


		//InitMyStr(nxp1,ny,&Hx);
		//InitMyStr(nx,nyp1,&Hy);
		//InitMyStr(nxp1,nyp1,&Ez);

		/*---------------< X >-------------------*/
		InitMyStr(nxp1,ny,&Hx_s);
		//InitMyStr(nxp1,ny,&Hx_i);
		//
		//InitMyStr(nxp1,ny,&Hx_i_pre);
		/*----------------------------------*/

		/*---------------< Y >-----------------*/
		InitMyStr(nx,nyp1,&Hy_s);
		
		//InitMyStr(nx,nyp1,&Hy_i);
		//InitMyStr(nx,nyp1,&Hy_i_pre);
		/*----------------------------------*/


		/*--------------< Z >-----------------*/
		InitMyStr(nxp1,nyp1,&Ez_s);
		//InitMyStr(nxp1,nyp1,&Ez_i);

		//InitMyStr(nxp1,nyp1,&Ez_i_pre);

		InitMyStr(nxp1,nyp1,&Vez);
		InitMyStr(nxp1,nyp1,&Ceez);
		InitMyStr(nxp1,nyp1,&Cevz);
		InitMyStr(nxp1,nyp1,&Cehz);
		/*----------------------------------*/

	}
	if(IfWithDensity){
		InitMyStr(m*nxp1,m*nyp1,&Erms);
		InitMyStr(m*nxp1,m*nyp1,&ne);
		InitMyStr(m*nxp1,m*nyp1,&ne_pre);
		InitMyStr(m*nxp1,m*nyp1,&beta);
		if(niutype==5){
			InitMyStr(m*nxp1,m*nyp1,&Nu_c);
			ResetStructDataR(Nu_c,5e9*760.0);
			if(IsTEx){
				InitMyStr(nxp1,nyp1,&Cvvz);
				InitMyStr(nxp1,nyp1,&Cvez);
			}
			if(IsTMx){
				InitMyStr(nx,nyp1,&Cvvx);
				InitMyStr(nx,nyp1,&Cvex);
				InitMyStr(nxp1,ny,&Cvvy);
				InitMyStr(nxp1,ny,&Cvey);
			}
			//InitMyStr(m*nxp1,m*nyp1,&Nu_c_pre);
		}
	}
}

#endif

