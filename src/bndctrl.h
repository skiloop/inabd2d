#ifndef BOUNDARY_CONTROL_INCLUDE
#define BOUNDARY_CONTROL_INCLUDE

int IsPMLxn = 0;
int IsPMLxp = 0;
int IsPMLyn = 0;
int IsPMLyp = 0;

int IsAnySidePML = 0;

int PMLCellxn = 0;
int PMLCellxp = 0;
int PMLCellyn = 0;
int PMLCellyp = 0;

int pml_order = 2;
MyDataF pml_R_0 = 1e-8;

int pis,pie;
int pjs,pje;
void InitBndCtrl(){
	IsAnySidePML = 1;
	PMLCellxp=nbound;
	PMLCellyp=nbound;
	PMLCellxn=nbound;
	PMLCellyn=nbound;
	IsPMLxp=1;
	IsPMLyp=1;
	IsPMLxn=1;
	IsPMLyn=1;
	
	pis = PMLCellxn;
	pie = nx-PMLCellxp;

	pjs = PMLCellyn;
	pje = ny-PMLCellyp;
}

#endif

