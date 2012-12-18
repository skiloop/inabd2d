void AddBndMurXN(MyStruct nf,MyStruct pf,int nbd){
	int i=0,j=0;
	MyDataF C=(c*dt-dx)/(c*dt+dx);
	for(i=nbd-1;i>=0;i--)
		for(j=nbd;j<nf.ny-nbd;j++)
			nf.data[i*nf.ny+j]=pf.data[(i+1)*pf.ny+j]+C*(nf.data[(i+1)*nf.ny+j]-nf.data[i*nf.ny+j]);
}
void AddBndMurXP(MyStruct nf,MyStruct pf,int nbd){
	int i=0,j=0;
	MyDataF C=(c*dt-dx)/(c*dt+dx);
	for(i=nf.nx-nbd;i<nf.nx;i++)
		for(j=nbd;j<nf.ny-nbd;j++)
			nf.data[i*nf.ny+j]=pf.data[(i-1)*pf.ny+j]+C*(nf.data[(i-1)*nf.ny+j]-nf.data[i*nf.ny+j]);
}
void AddBndMurYN(MyStruct nf,MyStruct pf,int nbd){
	int i=0,j=0;
	MyDataF C=(c*dt-dy)/(c*dt+dy);
	for(i=nbd;i<nf.nx-nbd;i++)
		for(j=nbd-1;j>=0;j--)
			nf.data[i*nf.ny+j]=pf.data[i*pf.ny+j+1]+C*(nf.data[i*nf.ny+j+1]-nf.data[i*nf.ny+j]);
}
void AddBndMurYP(MyStruct nf,MyStruct pf,int nbd){
	int i=0,j=0;
	MyDataF C=(c*dt-dy)/(c*dt+dy);
	for(i=nbd;i<nf.nx-nbd;i++)
		for(j=nf.ny-nbd;j<nf.ny;j++)
			nf.data[i*nf.ny+j]=pf.data[i*pf.ny+j-1]+C*(nf.data[i*nf.ny+j-1]-nf.data[i*nf.ny+j]);
}
void AddBndMurXNYN(MyStruct nf,MyStruct pf,int nbd){
	int i=0,j=0;
	MyDataF C=(c*dt-1.414*dy)/(c*dt+1.414*dy);
	for(i=nbd-1;i>=0;i--)
		for(j=nbd-1;j>=0;j--)
			nf.data[i*nf.ny+j]=pf.data[(i+1)*pf.ny+j+1]+C*(nf.data[(i+1)*nf.ny+j+1]-nf.data[i*nf.ny+j]);
}
void AddBndMurXPYN(MyStruct nf,MyStruct pf,int nbd){
	int i=0,j=0;
	MyDataF C=(c*dt-1.414*dy)/(c*dt+1.414*dy);
	for(i=nf.nx-nbd;i<nf.nx;i++)
		for(j=nbd-1;j>=0;j--)
			nf.data[i*nf.ny+j]=pf.data[(i-1)*pf.ny+j+1]+C*(nf.data[(i-1)*nf.ny+j+1]-nf.data[i*nf.ny+j]);
}
void AddBndMurXPYP(MyStruct nf,MyStruct pf,int nbd){
	int i=0,j=0;
	MyDataF C=(c*dt-1.414*dy)/(c*dt+1.414*dy);
	for(i=nf.nx-nbd;i<nf.nx;i++)
		for(j=nf.ny-nbd;j<nf.ny;j++)
			nf.data[i*nf.ny+j]=pf.data[(i-1)*pf.ny+j-1]+C*(nf.data[(i-1)*nf.ny+j-1]-nf.data[i*nf.ny+j]);
}
void AddBndMurXNYP(MyStruct nf,MyStruct pf,int nbd){
	int i=0,j=0;
	MyDataF C=(c*dt-1.414*dy)/(c*dt+1.414*dy);
	for(i=nbd-1;i>=0;i--)
		for(j=nf.ny-nbd;j<nf.ny;j++)
			nf.data[i*nf.ny+j]=pf.data[(i+1)*pf.ny+j-1]+C*(nf.data[(i+1)*nf.ny+j-1]-nf.data[i*nf.ny+j]);
}
