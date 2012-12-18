

#include <stdio.h>
#include <stdlib.h>

void calErmsAtCoarsePoint()
{
	int i,j,sizex,sizey,hm;
	int index, ind;
	MDOUBLE tmp;
	hm = m/2;
	for(i=tpis;i<=tpie;i++){
		for(j=tpjs;j<=tpje;j++){
			ind = (i*m+hm)*Erms.ny+j*m;
			tmp = 0.25*(Ey_s.data[i*Ey_s.ny+j]+Ey_s.data[(i-1)*Ey_s.ny+j]+Ey_s.data[i*Ey_s.ny+j+1]+Ey_s.data[(i-1)*Ey_s.ny+j+1]);
			Erms.data[ind] += Ex_s.data[i*Ex_s.ny+j]*Ex_s.data[i*Ex_s.data[i*Ex_s.ny+j]+tmp*tmp;
			ind = (i*m)*Erms.ny+j*m+hm;
			tmp = 0.25*(Ex_s.data[i*Ex_s.ny+j]+Ex_s.data[(i+1)*Ex_s.ny+j]+Ex_s.data[i*Ex_s.ny+j-1]+Ex_s.data[(i+1)*Ex_s.ny+j-1]);
			Erms.data[ind] += Ey_s.data[i*Ey_s.ny+j]*Ey_s.data[i*Ey_s.data[i*Ey_s.ny+j]+tmp*tmp;
		}
	}
}

