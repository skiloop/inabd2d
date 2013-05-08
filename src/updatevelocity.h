/*
*		This file defines functions that update velcocity.
*
*
*
*
*/


void UpdateVelocity(){
	int i,j,index;
	//for format five and for
	if(niutype==5||niutype==4){
		if(IsTMx){
			for(i=0;i<Vex.nx;i++){
				for(j=0;j<Vex.ny;j++){
					index=i*Vex.ny+j;
					Vex.data[index]=Cvvx.data[index]*Vex.data[index]-Cvex.data[index]*(Ex_s.data[index]+Ex_s_pre.data[index]);
				}
			}
			for(i=0;i<Vey.nx;i++){
				for(j=0;j<Vey.ny;j++){
					index=i*Vey.ny+j;
					Vey.data[index]=Cvvy.data[index]*Vey.data[index]-Cvey.data[index]*(Ey_s.data[index]+Ey_s_pre.data[index]);
				}
			}
		}
		if(IsTEx){
			for(i=0;i<Vez.nx;i++){
				for(j=0;j<Vez.ny;j++){
					index = i*Vez.ny+j;
					Vez.data[index]=Cvvz.data[index]*Vez.data[index]-Cvez.data[index]*(Ez_s.data[index]+Ez_s_pre.data[index]);
				}
			}
		}
		
	}else{//for other format
		if(IsTMx){
			for(i=0;i<Vex.nx;i++){
				for(j=0;j<Vex.ny;j++){
					Vex.data[i*Vex.ny+j]=alpha*SD(Vex,i,j)-Cve*(SD(Ex_s,i,j)+SD(Ex_s_pre,i,j));
				}
			}
			for(i=0;i<Vey.nx;i++){
				for(j=0;j<Vey.ny;j++){
					Vey.data[i*Vey.ny+j]=alpha*SD(Vey,i,j)-Cve*(SD(Ey_s_pre,i,j)+SD(Ey_s,i,j));
				}
			}
		}
		if(IsTEx){
			for(i=0;i<Vez.nx;i++){
				for(j=0;j<Vez.ny;j++){
					index = i*Vez.ny+j;
					Vez.data[index]=alpha*Vez.data[index]-Cve*(Ez_s_pre.data[index]+Ez_s.data[index]);
				}
			}
		}
	}
} 
