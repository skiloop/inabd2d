#include"murbndtofields.h"
void MurBndToStru(MyStruct mystr,int nbd){
	// Append Mur absord boundary to Field

	AddBndMurXN(mystr,nbd);
	AddBndMurYN(mystr,nbd);
	AddBndMurXP(mystr,nbd);
	AddBndMurYP(mystr,nbd);
	AddBndMurXNYN(mystr,nbd);
	AddBndMurXNYP(mystr,nbd);
	AddBndMurXPYN(mystr,nbd);
	AddBndMurXPYP(mystr,nbd);
}


void AddBndEFieldTEx(){
	MurBndToStru(Ez_s,nbound);	
}
void AddBndEFieldTMx(){
	MurBndToStru(Ex_s,nbound);
	MurBndToStru(Ey_s,nbound);
}
void AddBndMFieldTEx(){
	MurBndToStru(Hx_s,nbound);
	MurBndToStru(Hy_s,nbound);
	
}
void AddBndMFieldTMx(){
	MurBndToStru(Hz_s,nbound);
	
}
void MurBoundEField(){
	if(IsTEx)
		AddBndEFieldTEx();
	if(IsTMx)
		AddBndEFieldTMx();
}
void MurBoundMField(){
	if(IsTEx)
		AddBndMFieldTEx();
	if(IsTMx)
		AddBndMFieldTMx();
}