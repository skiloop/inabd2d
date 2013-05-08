#ifndef SCATTER_TOTAL_BOUNDARY_INCLUDE
#define SCATTER_TOTAL_BOUNDARY_INCLUDE

#ifndef SCATTER_FIELD_DOMAIN_BND_SIZE
#define SCATTER_FIELD_DOMAIN_BND_SIZE 4
#endif

#include"bndctrl.h"
#include"initbnd.h"
#include"source.h"
void IntTtlFldDmnBnd(){
	//this function define where the total field domain is
	//Note: pis,pjs,pie,pje must be defined and initialzed 
	//		before calling this function
	tpis = pis + SCATTER_FIELD_DOMAIN_BND_SIZE;
	tpjs = pjs + SCATTER_FIELD_DOMAIN_BND_SIZE;
	tpie = pie - SCATTER_FIELD_DOMAIN_BND_SIZE;
	tpje = pje - SCATTER_FIELD_DOMAIN_BND_SIZE;
	printf("\n*************************************************************\n");
	printf("\ntpis\t=\t%d\ntpjs\t=\t%d\ntpie\t=\t%d\ntpje\t=\t%d\n",tpis,tpjs,tpie,tpje);
	CalDelays();
}






void Adjust_E_Field(MyDataF pre_t)
{//adjust TEx electric field which at total fields boundary
	int ind;
	int index;
	//bound yn,yp
	
	for(ind = tpis;ind<=tpie;ind++){
		index = ind*Ez_s.ny;
		//bottom  Cehz.data[index+tpjs]
		Ez_s.data[index+tpjs]+=0;//Cehz.data[index+tpjs]/dy*Hx0*Source(pre_t-Dhxb[ind-tpis]*t_per_cell);
		//top
		Ez_s.data[index+tpje]-=0;//Cehz.data[index+tpje]/dy*Hx0*Source(pre_t-Dhxt[ind-tpis]*t_per_cell);
	}
	//bound xn,xp Cehz.data[index+tpje]
	for(ind = tpjs;ind<=tpje;ind++){
		index = tpis*Ez_s.ny+ind;
		//left sideCehz.data[index]
		Ez_s.data[index]-=Cehz.data[index]/dx*Hy0*Source(pre_t-Dhyl[ind-tpjs]*t_per_cell);
		index = tpie*Ez_s.ny+ind;
		//right
		Ez_s.data[index]+=Cehz.data[index]/dx*Hy0*Source(pre_t-Dhyr[ind-tpjs]*t_per_cell);
	}

}





//
void Adjust_M_Field(MyDataF pre_t){
	//adjust magnetic field of TMx in total field boundary
	int ind;
	MyDataF tmp;
	//Adjust Hy
	for(ind=tpjs;ind<=tpje;ind++){
		//left
		Hy_s.data[(tpis-1)*Hy_s.ny+ind]	-=	Chyez*Ez0*Source(pre_t-Dezl[ind-tpjs]*t_per_cell);
		//right
		Hy_s.data[tpie*Hy_s.ny+ind]		+=	Chyez*Ez0*Source(pre_t-Dezr[ind-tpjs]*t_per_cell);
	}
	//Adjust Hy
	for(ind=tpis;ind<=tpie;ind++){
		MyDataF Einc;
		
		//bottom
		Einc=Ez0*Source(pre_t-Dezb[ind-tpis]*t_per_cell);
		tmp = Chxez*Einc;
		Hx_s.data[ind*Hx_s.ny+tpjs-1]	-=	0;//tmp;//0;
		//top
		Hx_s.data[ind*Hx_s.ny+tpje]		+=	0;//Chxez*Ez0*Source(pre_t-Dezt[ind-tpis]*t_per_cell);
	}
}


#endif
