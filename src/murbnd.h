#include"murbndtofields.h"
#include"murnontol.h"
void TExMFieldMur(){
	//Hx_s
	AddBndMurXN(Hx_s,Hx_s_pre,nbound);
	AddBndMurXP(Hx_s,Hx_s_pre,nbound);
	AddBndMurYN(Hx_s,Hx_s_pre,nbound);
	AddBndMurYP(Hx_s,Hx_s_pre,nbound);

	AddBndMurXNYN(Hx_s,Hx_s_pre,nbound);
	AddBndMurXNYP(Hx_s,Hx_s_pre,nbound);
	AddBndMurXPYN(Hx_s,Hx_s_pre,nbound);
	AddBndMurXPYP(Hx_s,Hx_s_pre,nbound);

	//Hy_s
	AddBndMurXN(Hy_s,Hy_s_pre,nbound);
	AddBndMurXP(Hy_s,Hy_s_pre,nbound);
	AddBndMurYN(Hy_s,Hy_s_pre,nbound);
	AddBndMurYP(Hy_s,Hy_s_pre,nbound);

	AddBndMurXNYN(Hy_s,Hy_s_pre,nbound);
	AddBndMurXNYP(Hy_s,Hy_s_pre,nbound);
	AddBndMurXPYN(Hy_s,Hy_s_pre,nbound);
	AddBndMurXPYP(Hy_s,Hy_s_pre,nbound);
}
void TExEFieldMur(){
	AddBndMurXN(Ez_s,Ez_s_pre,nbound);
	AddBndMurXP(Ez_s,Ez_s_pre,nbound);
	AddBndMurYN(Ez_s,Ez_s_pre,nbound);
	AddBndMurYP(Ez_s,Ez_s_pre,nbound);

	AddBndMurXNYN(Ez_s,Ez_s_pre,nbound);
	AddBndMurXNYP(Ez_s,Ez_s_pre,nbound);
	AddBndMurXPYN(Ez_s,Ez_s_pre,nbound);
	AddBndMurXPYP(Ez_s,Ez_s_pre,nbound);
}
void TMxEFieldMur(){
	//Ex_s
	AddBndMurXN(Ex_s,Ex_s_pre,nbound);
	AddBndMurXP(Ex_s,Ex_s_pre,nbound);
	AddBndMurYN(Ex_s,Ex_s_pre,nbound);
	AddBndMurYP(Ex_s,Ex_s_pre,nbound);

	AddBndMurXNYN(Ex_s,Ex_s_pre,nbound);
	AddBndMurXNYP(Ex_s,Ex_s_pre,nbound);
	AddBndMurXPYN(Ex_s,Ex_s_pre,nbound);
	AddBndMurXPYP(Ex_s,Ex_s_pre,nbound);

	//Ey_s
	AddBndMurXN(Ey_s,Ey_s_pre,nbound);
	AddBndMurXP(Ey_s,Ey_s_pre,nbound);
	AddBndMurYN(Ey_s,Ey_s_pre,nbound);
	AddBndMurYP(Ey_s,Ey_s_pre,nbound);

	AddBndMurXNYN(Ey_s,Ey_s_pre,nbound);
	AddBndMurXNYP(Ey_s,Ey_s_pre,nbound);
	AddBndMurXPYN(Ey_s,Ey_s_pre,nbound);
	AddBndMurXPYP(Ey_s,Ey_s_pre,nbound);
}
void TMxMFieldMur(){
	AddBndMurXN(Hz_s,Hz_s_pre,nbound);
	AddBndMurXP(Hz_s,Hz_s_pre,nbound);
	AddBndMurYN(Hz_s,Hz_s_pre,nbound);
	AddBndMurYP(Hz_s,Hz_s_pre,nbound);

	AddBndMurXNYN(Hz_s,Hz_s_pre,nbound);
	AddBndMurXNYP(Hz_s,Hz_s_pre,nbound);
	AddBndMurXPYN(Hz_s,Hz_s_pre,nbound);
	AddBndMurXPYP(Hz_s,Hz_s_pre,nbound);}

void MFieldMurAB(){
	if(IsTEx)
		TExMFieldMur();
	if(IsTMx)
		TMxMFieldMur();
}
void EFieldMurAB(){
	if(IsTEx)
		TExEFieldMur();
	if(IsTMx)
		TMxEFieldMur();
}



