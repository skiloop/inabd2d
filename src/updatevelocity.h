/*
*		This file defines functions that update velcocity.
*
*
*
*
*/


void UpdateVelocity(){
	int i,j,index;
	if(IsTMx){
		for(i=0;i<Vex.nx;i++)
			for(j=0;j<Vex.ny;j++){
				Vex.data[i*Vex.ny+j]=alpha*SD(Vex,i,j)-Cve*(SD(Ex_s,i,j)+SD(Ex_s_pre,i,j));
			}
		for(i=0;i<Vey.nx;i++)
			for(j=0;j<Vey.ny;j++){
				Vey.data[i*Vey.ny+j]=alpha*SD(Vey,i,j)-Cve*(SD(Ey_s_pre,i,j)+SD(Ey_s,i,j));
			}
	}
	if(IsTEx){
		for(i=0;i<Vez.nx;i++)
			for(j=0;j<Vez.ny;j++){
				index = i*Vez.ny+j;
				Vez.data[index]=alpha*Vez.data[index]-Cve*(Ez_s_pre.data[index]+Ez_s.data[index]);

				/*if(Vez.data[i*Vez.ny+j]>1e7||Vez.data[i*Vez.ny+j]<-1e7){
					printf("Vez[%d,%d]=%f",i,j,Vez.data[i*Vez.ny+j]);
					//system("pause");
				}*/

			}
	}
} 
