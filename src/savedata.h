
#define GetD(k,i,j) k.data[(i)*k.ny+j]


/******************************************************/

void SaveData(MyStruct data,char* filename,FILE *hfile);



/******************************************************/
void SaveToFile(){

	FILE *hfile;
	char dlist[15]="DataList.dat";

	if((hfile=fopen(dlist,"w"))==NULL){
		printf("Cannot open file: \t%s\n",dlist);
		exit(0);
	}
	fprintf(hfile,"DataName\tnx\tny\n");
	
	if(IsTMx){
		//SaveData(Ex,"Ex.dat",hfile);
		//SaveData(Ey,"Ey.dat",hfile);
		//SaveData(Hz,"Hz.dat",hfile);

		SaveData(Ex_s,"Ex_s.dat",hfile);
		SaveData(Ey_s,"Ey_s.dat",hfile);
		SaveData(Hz_s,"Hz_s.dat",hfile);

		//SaveData(Ex_i,"Ex_i.dat",hfile);
		//SaveData(Ey_i,"Ey_i.dat",hfile);
		//SaveData(Hz_i,"Hz_i.dat",hfile);

		SaveData(Vex,"Vex.dat",hfile);
		SaveData(Vey,"Vey.dat",hfile);
		
	}
	if(IsTEx){
		//SaveData(Hx,"Hx.dat",hfile);
		//SaveData(Hy,"Hy.dat",hfile);
		//SaveData(Ez,"Ez.dat",hfile);

		SaveData(Hx_s,"Hx_s.dat",hfile);
		SaveData(Hy_s,"Hy_s.dat",hfile);
		SaveData(Ez_s,"Ez_s.dat",hfile);
/*
		SaveData(Hx_i,"Hx_i.dat",hfile);
		SaveData(Hy_i,"Hy_i.dat",hfile);
		SaveData(Ez_i,"Ez_i.dat",hfile);
*/
		SaveData(Vez,"Vez.dat",hfile);
		
	}
	SaveData(ne,"ne.dat",hfile);
	fclose(hfile);
}

void SaveData(MyStruct data,char* filename,FILE *hfile){
	int i,j;
	FILE *fp;

	if((fp=fopen(filename,"w"))==NULL){
			printf("Cannot create file: \t%s\n",filename);
			exit(0);
		}
	
	fprintf(hfile,"%s\t%d\t%d\n",filename,data.nx,data.ny);
	
	printf("\n%s\n",filename);
	//PrintData(data);
	//fwrite(data.data,sizeof(MyDataF),data.nx*data.ny,fp);
	for(j=0;j<data.ny;j++){
		for(i=0;i<data.nx;i++)
			fprintf(fp,"%2.3E ",GetD(data,i,j));
		fprintf(fp,"\n");
	}
	
	fclose(fp);
}
