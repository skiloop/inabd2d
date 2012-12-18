
void TExMFieldMur2(){
	//Hx
	AddBndMurXN(Hx,Hx_pre,nbound);
	AddBndMurXP(Hx,Hx_pre,nbound);
	AddBndMurYN(Hx,Hx_pre,nbound);
	AddBndMurYP(Hx,Hx_pre,nbound);

	AddBndMurXNYN(Hx,Hx_pre,nbound);
	AddBndMurXNYP(Hx,Hx_pre,nbound);
	AddBndMurXPYN(Hx,Hx_pre,nbound);
	AddBndMurXPYP(Hx,Hx_pre,nbound);

	//Hy
	AddBndMurXN(Hy,Hy_pre,nbound);
	AddBndMurXP(Hy,Hy_pre,nbound);
	AddBndMurYN(Hy,Hy_pre,nbound);
	AddBndMurYP(Hy,Hy_pre,nbound);

	AddBndMurXNYN(Hy,Hy_pre,nbound);
	AddBndMurXNYP(Hy,Hy_pre,nbound);
	AddBndMurXPYN(Hy,Hy_pre,nbound);
	AddBndMurXPYP(Hy,Hy_pre,nbound);
}
void TExEFieldMur2(){
	AddBndMurXN(Ez,Ez_pre,nbound);
	AddBndMurXP(Ez,Ez_pre,nbound);
	AddBndMurYN(Ez,Ez_pre,nbound);
	AddBndMurYP(Ez,Ez_pre,nbound);

	AddBndMurXNYN(Ez,Ez_pre,nbound);
	AddBndMurXNYP(Ez,Ez_pre,nbound);
	AddBndMurXPYN(Ez,Ez_pre,nbound);
	AddBndMurXPYP(Ez,Ez_pre,nbound);
}
void TMxEFieldMur2(){
	//Ex
	AddBndMurXN(Ex,Ex_pre,nbound);
	AddBndMurXP(Ex,Ex_pre,nbound);
	AddBndMurYN(Ex,Ex_pre,nbound);
	AddBndMurYP(Ex,Ex_pre,nbound);

	AddBndMurXNYN(Ex,Ex_pre,nbound);
	AddBndMurXNYP(Ex,Ex_pre,nbound);
	AddBndMurXPYN(Ex,Ex_pre,nbound);
	AddBndMurXPYP(Ex,Ex_pre,nbound);

	//Ey
	AddBndMurXN(Ey,Ey_pre,nbound);
	AddBndMurXP(Ey,Ey_pre,nbound);
	AddBndMurYN(Ey,Ey_pre,nbound);
	AddBndMurYP(Ey,Ey_pre,nbound);

	AddBndMurXNYN(Ey,Ey_pre,nbound);
	AddBndMurXNYP(Ey,Ey_pre,nbound);
	AddBndMurXPYN(Ey,Ey_pre,nbound);
	AddBndMurXPYP(Ey,Ey_pre,nbound);
}
void TMxMFieldMur2(){
	AddBndMurXN(Hz,Hz_pre,nbound);
	AddBndMurXP(Hz,Hz_pre,nbound);
	AddBndMurYN(Hz,Hz_pre,nbound);
	AddBndMurYP(Hz,Hz_pre,nbound);

	AddBndMurXNYN(Hz,Hz_pre,nbound);
	AddBndMurXNYP(Hz,Hz_pre,nbound);
	AddBndMurXPYN(Hz,Hz_pre,nbound);
	AddBndMurXPYP(Hz,Hz_pre,nbound);
}

void MFieldMurAB2(){
	if(IsTEx)
		TExMFieldMur2();
	if(IsTMx)
		TMxMFieldMur2();
}
void EFieldMurAB2(){
	if(IsTEx)
		TExEFieldMur2();
	if(IsTMx)
		TMxEFieldMur2();
}



