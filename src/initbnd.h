#ifndef _INIT_BND_H_
#define _INIT_BND_H_


void InitPMLBoundTEx(){
	int i,j,ind;
	
	MDOUBLE	*rho_e;
	MDOUBLE	*rho_m;
	
	COORDINATES	sigma_pey_yn;
	COORDINATES	sigma_pmy_yn;

	COORDINATES	sigma_pex_xp;
	COORDINATES	sigma_pmx_xp;

	COORDINATES	sigma_pex_xn;
	COORDINATES	sigma_pmx_xn;

	COORDINATES	sigma_pey_yp;
	COORDINATES	sigma_pmy_yp;

	MDOUBLE	sigma_max;

	printf("Initializing pml boundary conditions 2d TEx...\n");
	//printf("nxp1=%d\nnym1=%d\nny=%d\n",nxp1,nym1,ny);
	init_coordinates(&Ezx_xn,PMLCellxn,nym1,0.0);
	init_coordinates(&Ezx_xp,PMLCellxp,nym1,0.0);
	init_coordinates(&Ezx_yn,nxm1-PMLCellxn-PMLCellxp,PMLCellyn,0.0);
	init_coordinates(&Ezx_yp,nxm1-PMLCellxn-PMLCellxp,PMLCellyp,0.0);	
	init_coordinates(&Ezy_xn,PMLCellxn,nym1-PMLCellyn-PMLCellyp,0.0);
	init_coordinates(&Ezy_xp,PMLCellxp,nym1-PMLCellyn-PMLCellyp,0.0);
	init_coordinates(&Ezy_yn,nxm1,PMLCellyn,0.0);
	init_coordinates(&Ezy_yp,nxm1,PMLCellyp,0.0);	

	if(IsPMLxn){

		rho_e = (MDOUBLE*)malloc(PMLCellxn*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellxn*sizeof(MDOUBLE));

		init_coordinates(&sigma_pex_xn,PMLCellxn,nym1,0.0);
		init_coordinates(&sigma_pmx_xn,PMLCellxn,nym1,0.0);		
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxn);

		for(i=PMLCellxn;i>0;--i){
			rho_e[PMLCellxn-i] = (i-0.75)/PMLCellxn;
			rho_m[PMLCellxn-i] = (i-0.25)/PMLCellxn;
		}

		for(i=0;i<PMLCellxn;i++)
			for(j=0;j<nym1;++j){
				sigma_pex_xn.data[i*nym1+j] = sigma_max*pow(rho_e[i],pml_order);
				sigma_pmx_xn.data[i*nym1+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
			}


		init_coordinates(&Chyh_xn,PMLCellxn,nym1,0.0);
		init_coordinates(&Chyez_xn,PMLCellxn,nym1,0.0);

		init_coordinates(&Cezxe_xn,Ezx_xn.nx,Ezx_xn.ny,0.0);
		init_coordinates(&Cezxhy_xn,Ezx_xn.nx,Ezx_xn.ny,0.0);

		init_coordinates(&Cezye_xn ,Ezy_xn.nx,Ezy_xn.ny,1.0);
		init_coordinates(&Cezyhx_xn,Ezy_xn.nx,Ezy_xn.ny,-dt/(dy*eps_0));
		
		ind = 0;
		for(i=0;i<PMLCellxn;++i)
			for(j=0;j<nym1;j++){
				/* Coefficients updating Hx */
				Chyh_xn.data[ind] = (2*mu_0-dt*sigma_pmx_xn.data[ind])/(2*mu_0+dt*sigma_pmx_xn.data[ind]);
				Chyez_xn.data[ind] = (2*dt/dx)/(2*mu_0+dt*sigma_pmx_xn.data[ind]);
				
				/* Coefficients updating Ezx */
				Cezxe_xn.data[ind] = (2*eps_0-dt*sigma_pex_xn.data[ind])/(2*eps_0+dt*sigma_pex_xn.data[ind]);
				Cezxhy_xn.data[ind] = (2*dt/dx)/(2*eps_0+dt*sigma_pex_xn.data[ind]);

				ind++;
			}
		
		free(sigma_pmx_xn.data);
		free(sigma_pex_xn.data);

		free(rho_e);
		free(rho_m);
	}

	if(IsPMLyn){

		rho_e = (MDOUBLE*)malloc(PMLCellyn*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellyn*sizeof(MDOUBLE));

		init_coordinates(&sigma_pey_yn,nxm1,PMLCellyn,0.0);
		init_coordinates(&sigma_pmy_yn,nxm1,PMLCellyn,0.0);		
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyn);

		for(i=PMLCellyn;i>0;--i){
			rho_e[PMLCellyn-i] = (i-0.75)/PMLCellyn;
			rho_m[PMLCellyn-i] = (i-0.25)/PMLCellyn;
		}

		for(i=0;i<PMLCellyn;i++)
			for(j=0;j<nxm1;++j){
				sigma_pey_yn.data[j*PMLCellyn+i] = sigma_max*pow(rho_e[i],pml_order);
				sigma_pmy_yn.data[j*PMLCellyn+i] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
			}

		init_coordinates(&Chxh_yn,nxm1,PMLCellyn,0.0);
		init_coordinates(&Chxez_yn,nxm1,PMLCellyn,0.0);

		
		init_coordinates(&Cezye_yn,Ezy_yn.nx,Ezy_yn.ny,0.0);
		init_coordinates(&Cezyhx_yn,Ezy_yn.nx,Ezy_yn.ny,0.0);
		
		ind = 0;
		for(i=0;i<nxm1;++i)
			for(j=0;j<PMLCellyn;j++){
				/* Coefficients updating Hx */
				Chxh_yn.data[ind] = (2*mu_0-dt*sigma_pmy_yn.data[ind])/(2*mu_0+dt*sigma_pmy_yn.data[ind]);
				Chxez_yn.data[ind] = -(2*dt/dy)/(2*mu_0+dt*sigma_pmy_yn.data[ind]);
			
				/* Coefficients updating Ezy */
				Cezye_yn.data[ind] = (2*eps_0-dt*sigma_pey_yn.data[ind])/(2*eps_0+dt*sigma_pey_yn.data[ind]);
				Cezyhx_yn.data[ind] = -(2*dt/dy)/(2*eps_0+dt*sigma_pey_yn.data[ind]);
				ind++;
			}

		
		init_coordinates(&Cezxe_yn,Ezx_yn.nx,Ezx_yn.ny,1.0);
		init_coordinates(&Cezxhy_yn,Ezx_yn.nx,Ezx_yn.ny,dt/(dx*eps_0));
		
		free(sigma_pmy_yn.data);
		free(sigma_pey_yn.data);

		free(rho_e);
		free(rho_m);
	}

	if(IsPMLxp){

		rho_e = (MDOUBLE*)malloc(PMLCellxp*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellxp*sizeof(MDOUBLE));

		init_coordinates(&sigma_pex_xp,PMLCellxp,nym1,0.0);
		init_coordinates(&sigma_pmx_xp,PMLCellxp,nym1,0.0);
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxp);

		for(i=0;i<PMLCellyp;++i){
			rho_e[i] = (i+1-0.75)/PMLCellxp;
			rho_m[i] = (i+1-0.25)/PMLCellxp;
		}

		for(i=0;i<nym1;i++)
			for(j=0;j<PMLCellxp;++j){
				sigma_pex_xp.data[j*nym1+i] = sigma_max*pow(rho_e[j],pml_order);
				sigma_pmx_xp.data[j*nym1+i] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
			}
		
		init_coordinates(&Chyh_xp,PMLCellxp,nym1,0.0);
		init_coordinates(&Chyez_xp,PMLCellxp,nym1,0.0);

		init_coordinates(&Cezxe_xp,Ezx_xp.nx,Ezx_xp.ny,0.0);
		init_coordinates(&Cezxhy_xp,Ezx_xp.nx,Ezx_xp.ny,0.0);		

		ind = 0;
		for(i=0;i<PMLCellxp;++i)
			for(j=0;j<nym1;j++){
				/* Coefficients updating Hy */
				Chyh_xp.data[ind] = (2*mu_0-dt*sigma_pmx_xp.data[ind])/(2*mu_0+dt*sigma_pmx_xp.data[ind]);
				Chyez_xp.data[ind] = (2*dt/dx)/(2*mu_0+dt*sigma_pmx_xp.data[ind]);
				
				/* Coefficients updating Ezx */
				Cezxe_xp.data[ind] = (2*eps_0-dt*sigma_pex_xp.data[ind])/(2*eps_0+dt*sigma_pex_xp.data[ind]);
				Cezxhy_xp.data[ind] = (2*dt/dx)/(2*eps_0+dt*sigma_pex_xp.data[ind]);

				ind++;
			}

		
		init_coordinates(&Cezye_xp,Ezy_xp.nx,Ezy_xp.ny,1.0);
		init_coordinates(&Cezyhx_xp,Ezy_xp.nx,Ezy_xp.ny,-dt/(dy*eps_0));
		
		free(sigma_pmx_xp.data);
		free(sigma_pex_xp.data);

		free(rho_e);
		free(rho_m);
	}
	if(IsPMLyp){

		rho_e = (MDOUBLE*)malloc(PMLCellyp*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellyp*sizeof(MDOUBLE));

		init_coordinates(&sigma_pey_yp,nxm1,PMLCellyp,0.0);
		init_coordinates(&sigma_pmy_yp,nxm1,PMLCellyp,0.0);
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyp);

		for(i=0;i<PMLCellyp;++i){
			rho_e[i] = (i+1-0.75)/PMLCellyp;
			rho_m[i] = (i+1-0.25)/PMLCellyp;
		}

		for(i=0;i<nxm1;i++)
			for(j=0;j<PMLCellyp;++j){
				sigma_pey_yp.data[i*PMLCellyp+j] = sigma_max*pow(rho_e[j],pml_order);
				sigma_pmy_yp.data[i*PMLCellyp+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
			}
		
		init_coordinates(&Chxh_yp,nxm1,PMLCellyp,0.0);
		init_coordinates(&Chxez_yp,nxm1,PMLCellyp,0.0);
		
		init_coordinates(&Cezye_yp, Ezy_yp.nx,Ezy_yp.ny,0.0);
		init_coordinates(&Cezyhx_yp,Ezy_yp.nx,Ezy_yp.ny,0.0);

		ind = 0;
		for(i=0;i<nxm1;++i)
			for(j=0;j<PMLCellyp;j++){
				/* Coefficients updating Hy */
				Chxh_yp.data[ind] = (2*mu_0-dt*sigma_pmy_yp.data[ind])/(2*mu_0+dt*sigma_pmy_yp.data[ind]);
				Chxez_yp.data[ind] = -(2*dt/dx)/(2*mu_0+dt*sigma_pmy_yp.data[ind]);
								
				/* Coefficients updating Ezy */
				Cezye_yp.data[ind] = (2*eps_0-dt*sigma_pey_yp.data[ind])/(2*eps_0+dt*sigma_pey_yp.data[ind]);
				Cezyhx_yp.data[ind] = -(2*dt/dy)/(2*eps_0+dt*sigma_pey_yp.data[ind]);
				ind++;
			}
		
		init_coordinates(&Cezxe_yp, Ezx_yp.nx,Ezx_yp.ny,1.0);
		init_coordinates(&Cezxhy_yp,Ezx_yp.nx,Ezx_yp.ny,dt/(dx*eps_0));
		
		free(sigma_pmy_yp.data);
		free(sigma_pey_yp.data);

		free(rho_e);
		free(rho_m);
	}
	

	printf("End of initializing pml boundary conditions 2d TEx...\n");
}
//void InitPMLBoundTMx(){
//
//	int i,j,ind;
//	MDOUBLE pmx,pex;
//	MDOUBLE	rho_e;
//	MDOUBLE	rho_m;
//
//	//COORDINATES	sigma_pex_xn;
//	//COORDINATES	sigma_pmx_xn;
//
//	//COORDINATES	sigma_pey_yn;
//	//COORDINATES	sigma_pmy_yn;
//
//	//COORDINATES	sigma_pex_xp;
//	//COORDINATES	sigma_pmx_xp;
//
//	//COORDINATES	sigma_pey_yp;
//	//COORDINATES	sigma_pmy_yp;
//	
//
//	MDOUBLE	sigma_max;
//
//	printf("Initializing pml boundary conditions 2d TMx...\n");
//	
//	init_coordinates(&Hzx_xn,PMLCellxn,ny,0.0);
//	init_coordinates(&Hzx_xp,PMLCellxp,ny,0.0);
//
//	init_coordinates(&Hzx_yn,nx-PMLCellxn-PMLCellxp,PMLCellyn,0.0);
//	init_coordinates(&Hzx_yp,nx-PMLCellxp-PMLCellxn,PMLCellyp,0.0);
//
//	init_coordinates(&Hzy_xn,PMLCellxn,ny-PMLCellyn-PMLCellyp,0.0);
//	init_coordinates(&Hzy_xp,PMLCellxp,ny-PMLCellyp-PMLCellyn,0.0);
//
//	init_coordinates(&Hzy_yn,nx,PMLCellyn,0.0);
//	init_coordinates(&Hzy_yp,nx,PMLCellyp,0.0);
//
///*=====================================================================================================================*/
//	if(IsPMLxn){
//		rho_e = (MDOUBLE*)malloc(PMLCellxn*sizeof(MDOUBLE));
//		rho_m = (MDOUBLE*)malloc(PMLCellxn*sizeof(MDOUBLE));		
//
//		init_coordinates(&sigma_pex_xn,PMLCellxn,ny,0.0);		
//		init_coordinates(&sigma_pmx_xn,PMLCellxn,ny,0.0);
//		
//		init_coordinates(&Ceye_xn,PMLCellxn,ny,0.0);		
//		init_coordinates(&Ceyhz_xn,PMLCellxn,ny,0.0);
//
//		init_coordinates(&Chzxh_xn, Hzx_xn.nx,Hzx_xn.ny,0.0);		
//		init_coordinates(&Chzxey_xn,Hzx_xn.nx,Hzx_xn.ny,0.0);
//
//		init_coordinates(&Chzyh_xn, Hzy_xn.nx,Hzy_xn.ny,1.0);	
//		init_coordinates(&Chzyex_xn,Hzy_xn.nx,Hzy_xn.ny,dt/(dy*mu_0));
//
//		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxn);
//
//		//for(i=PMLCellxn;i>0;--i){
//		//	rho_e[PMLCellxn-i] = (i-0.75)/PMLCellxn;
//		//	rho_m[PMLCellxn-i] = (i-0.25)/PMLCellxn;
//		//}
//
//		//for(i=0;i<PMLCellxn;i++)
//		//	for(j=0;j<ny;++j){
//		//		sigma_pex_xn.data[i*ny+j] = sigma_max*pow(rho_e[i],pml_order);
//		//		sigma_pmx_xn.data[i*ny+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
//		//	}
//
//
//
//		for(i=0;i<PMLCellxn;++i){
//			rho_e = (i+0.25)/PMLCellxn;
//			rho_m = (i+0.75)/PMLCellxn;
//			pex = sigma_max*pow(rho_e[i],pml_order);
//			pmx = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
//			for(j=0;j<ny;j++){
//				/* Coefficients updating Ey */
//				Ceye_xn.data[ind] = (2*eps_0-dt*pex)/(2*eps_0+dt*pex);
//				Ceyhz_xn.data[ind] = -(2*dt/dx)/(2*eps_0+dt*pex);
//			
//				/* Coefficients updating Hzx */
//				Chzxh_xn.data[ind] = (2*mu_0-dt*pmx)/(2*mu_0+dt*pmx);
//				Chzxey_xn.data[ind] = -(2*dt/dx)/(2*mu_0+dt*pmx);
//			}		
//		}
//		free(rho_e);
//		free(rho_m);
//	}
//
///*====================================================================================================================*/
//	if(IsPMLxp){
//
//		rho_e = (MDOUBLE*)malloc(PMLCellxp*sizeof(MDOUBLE));
//		rho_m = (MDOUBLE*)malloc(PMLCellxp*sizeof(MDOUBLE));
//
//		init_coordinates(&sigma_pex_xp,PMLCellxp,ny,0.0);		
//		init_coordinates(&sigma_pmx_xp,PMLCellxp,ny,0.0);
//		
//		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxp);
//
//		for(i=0;i<PMLCellxp;++i){
//			rho_e[i] = (i+1-0.75)/PMLCellxp;
//			rho_m[i] = (i+1-0.25)/PMLCellxp;
//		}
//
//		for(i=0;i<PMLCellxp;i++)
//			for(j=0;j<ny;++j){
//				sigma_pex_xp.data[i*ny+j] = sigma_max*pow(rho_e[i],pml_order);
//				sigma_pmx_xp.data[i*ny+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
//			}
//
//		init_coordinates(&Ceye_xp,PMLCellxp,ny,0.0);		
//		init_coordinates(&Ceyhz_xp,PMLCellxp,ny,0.0);
//
//		init_coordinates(&Chzxh_xp, Hzx_xp.nx,Hzx_xp.ny,0.0);		
//		init_coordinates(&Chzxey_xp,Hzx_xp.nx,Hzx_xp.ny,0.0);
//
//		init_coordinates(&Chzyh_xp, Hzy_xp.nx,Hzy_xp.ny,1.0);		
//		init_coordinates(&Chzyex_xp,Hzy_xp.nx,Hzy_xp.ny,dt/(dy*mu_0));	
//
//		ind = 0;
//		for(i=0;i<PMLCellxp;++i)
//			for(j=0;j<ny;j++){
//				/* Coefficients updating Ey */
//				Ceye_xp.data[ind] = (2*eps_0-dt*sigma_pex_xp.data[ind])/(2*eps_0+dt*sigma_pex_xp.data[ind]);
//				Ceyhz_xp.data[ind] = -(2*dt/dx)/(2*eps_0+dt*sigma_pex_xp.data[ind]);
//				
//				/* Coefficients updating Hzx */
//				Chzxh_xp.data[ind] = (2*mu_0-dt*sigma_pmx_xp.data[ind])/(2*mu_0+dt*sigma_pmx_xp.data[ind]);
//				Chzxey_xp.data[ind] = -(2*dt/dx)/(2*mu_0+dt*sigma_pmx_xp.data[ind]);
//
//				ind++;
//			}
//		
//		free(sigma_pmx_xp.data);
//		free(sigma_pex_xp.data);
//		
//		free(rho_e);
//		free(rho_m);
//		
//	}
///*=====================================================================================================================*/
//	if(IsPMLyn){
//		rho_e = (MDOUBLE*)malloc(PMLCellyn*sizeof(MDOUBLE));
//		rho_m = (MDOUBLE*)malloc(PMLCellyn*sizeof(MDOUBLE));
//
//		init_coordinates(&sigma_pey_yn,nx,PMLCellyn,0.0);		
//		init_coordinates(&sigma_pmy_yn,nx,PMLCellyn,0.0);
//		
//		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyn);
//
//		for(i=PMLCellyn;i>0;--i){
//			rho_e[PMLCellyn-i] = (i-0.75)/PMLCellyn;
//			rho_m[PMLCellyn-i] = (i-0.25)/PMLCellyn;
//		}
//
//		for(i=0;i<nx;i++)
//			for(j=0;j<PMLCellyn;++j){
//				sigma_pey_yn.data[i*PMLCellyn+j] = sigma_max*pow(rho_e[j],pml_order);
//				sigma_pmy_yn.data[i*PMLCellyn+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
//			}
//
//		init_coordinates(&Cexe_yn,nx,PMLCellyn,0.0);		
//		init_coordinates(&Cexhz_yn,nx,PMLCellyn,0.0);
//
//		init_coordinates(&Chzxh_yn, Hzx_yn.nx,Hzx_yn.ny,1.0);
//		init_coordinates(&Chzxey_yn,Hzx_yn.nx,Hzx_yn.ny,-dt/(dx*mu_0));
//
//		init_coordinates(&Chzyh_yn, Hzy_yn.nx,Hzy_yn.ny,0.0);		
//		init_coordinates(&Chzyex_yn,Hzy_yn.nx,Hzy_yn.ny,0.0);
//
//		ind = 0;
//		for(i=0;i<nx;++i)
//			for(j=0;j<PMLCellyn;j++){
//				/* Coefficients updating Ey */
//				Cexe_yn.data[ind] = (2*eps_0-dt*sigma_pey_yn.data[ind])/(2*eps_0+dt*sigma_pey_yn.data[ind]);
//				Cexhz_yn.data[ind] = (2*dt/dy)/(2*eps_0+dt*sigma_pey_yn.data[ind]);
//				
//				/* Coefficients updating Hzy */
//				Chzyh_yn.data[ind] = (2*mu_0-dt*sigma_pmy_yn.data[ind])/(2*mu_0+dt*sigma_pmy_yn.data[ind]);
//				Chzyex_yn.data[ind] = (2*dt/dy)/(2*mu_0+dt*sigma_pmy_yn.data[ind]);
//				ind++;
//			}
//		
//		free(sigma_pmy_yn.data);
//		free(sigma_pey_yn.data);
//
//		free(rho_e);
//		free(rho_m);
//	}
///*=====================================================================================================================*/
//	if(IsPMLyp){
//
//		rho_e = (MDOUBLE*)malloc(PMLCellyp*sizeof(MDOUBLE));
//		rho_m = (MDOUBLE*)malloc(PMLCellyp*sizeof(MDOUBLE));
//
//		init_coordinates(&sigma_pey_yp,nx,PMLCellyp,0.0);
//		init_coordinates(&sigma_pmy_yp,nx,PMLCellyp,0.0);
//		
//		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyp);
//
//		for(i=0;i<PMLCellyp;++i){
//			rho_e[i] = (i+1-0.75)/PMLCellyp;
//			rho_m[i] = (i+1-0.25)/PMLCellyp;
//		}
//
//		for(i=0;i<nx;i++)
//			for(j=0;j<PMLCellyp;++j){
//				sigma_pey_yp.data[i*PMLCellyp+j] = sigma_max*pow(rho_e[j],pml_order);
//				sigma_pmy_yp.data[i*PMLCellyp+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
//			}
//
//		init_coordinates(&Cexe_yp,nx,PMLCellyp,0.0);
//		init_coordinates(&Cexhz_yp,nx,PMLCellyp,0.0);
//
//		init_coordinates(&Chzxh_yp, Hzx_yp.nx,Hzx_yp.ny,1.0);
//		init_coordinates(&Chzxey_yp,Hzx_yp.nx,Hzx_yp.ny,-dt/(dx*mu_0));
//
//		init_coordinates(&Chzyh_yp,Hzy_yp.nx,Hzy_yp.ny,0.0);
//		init_coordinates(&Chzyex_yp,Hzy_yp.nx,Hzy_yp.ny,0.0);
//		ind = 0;
//		for(i=0;i<nx;++i)
//			for(j=0;j<PMLCellyp;j++){
//				/* Coefficients updating Ey */
//				Cexe_yp.data[ind] = (2*eps_0-dt*sigma_pey_yp.data[ind])/(2*eps_0+dt*sigma_pey_yp.data[ind]);
//				Cexhz_yp.data[ind] = (2*dt/dy)/(2*eps_0+dt*sigma_pey_yp.data[ind]);
//								
//				/* Coefficients updating Hzy */
//				Chzyh_yp.data[ind] = (2*mu_0-dt*sigma_pmy_yp.data[ind])/(2*mu_0+dt*sigma_pmy_yp.data[ind]);
//				Chzyex_yp.data[ind] = (2*dt/dy)/(2*mu_0+dt*sigma_pmy_yp.data[ind]);
//				ind++;
//			}
//	
//	free(sigma_pmy_yp.data);
//	free(sigma_pey_yp.data);
//
//	free(rho_e);
//	free(rho_m);
//	}
//	
//	printf("End of initializing pml boundary conditions 2d TMx...\n");
//}
//


void InitPMLBoundTMx(){

	int i,j,ind;
	MDOUBLE pmx,pex;
	MDOUBLE	rho_e;
	MDOUBLE	rho_m;
	MDOUBLE mu_d_eps = mu_0/eps_0;
	MDOUBLE	sigma_max;
	MDOUBLE tsigmax = -(pml_order+1)*eps_0*c*log(pml_R_0)/2;
	MDOUBLE dyepsdt = dy*eps_0/dt;
	MDOUBLE dxepsdt = dx*eps_0/dt;
	MDOUBLE dxmudt = dx*mu_0/dt;
	MDOUBLE dymudt = dy*mu_0/dt;
	MDOUBLE hdx = dx/2;
	MDOUBLE hdy = dy/2;
	printf("Initializing pml boundary conditions 2d TMx...\n");

	init_coordinates(&Hzx_xn,PMLCellxn,ny,0.0);
	init_coordinates(&Hzx_xp,PMLCellxp,ny,0.0);

	init_coordinates(&Hzx_yn,nx-PMLCellxn-PMLCellxp,PMLCellyn,0.0);
	init_coordinates(&Hzx_yp,nx-PMLCellxp-PMLCellxn,PMLCellyp,0.0);

	init_coordinates(&Hzy_xn,PMLCellxn,ny-PMLCellyn-PMLCellyp,0.0);
	init_coordinates(&Hzy_xp,PMLCellxp,ny-PMLCellyp-PMLCellyn,0.0);

	init_coordinates(&Hzy_yn,nx,PMLCellyn,0.0);
	init_coordinates(&Hzy_yp,nx,PMLCellyp,0.0);

	/*=====================================================================================================================*/
	if(IsPMLxn){

		init_coordinates(&Ceye_xn,PMLCellxn,ny,0.0);		
		init_coordinates(&Ceyhz_xn,PMLCellxn,ny,0.0);

		init_coordinates(&Chzxh_xn, Hzx_xn.nx,Hzx_xn.ny,0.0);		
		init_coordinates(&Chzxey_xn,Hzx_xn.nx,Hzx_xn.ny,0.0);

		init_coordinates(&Chzyh_xn, Hzy_xn.nx,Hzy_xn.ny,1.0);	
		init_coordinates(&Chzyex_xn,Hzy_xn.nx,Hzy_xn.ny,dt/(dy*mu_0));

		sigma_max = tsigmax/(dx*PMLCellxn);

		for(ind=i=0;i<PMLCellxn;++i){
			rho_e = 1-(i-0.75)/PMLCellxn;
			rho_m = 1-(i-0.25)/PMLCellxn;
			pex = sigma_max*pow(rho_e,pml_order);
			pmx = mu_d_eps*sigma_max*pow(rho_m,pml_order);
			for(j=0;j<ny;j++){
				/* Coefficients updating Ey */
				Ceye_xn.data[ind] = (2*eps_0-dt*pex)/(2*eps_0+dt*pex);
				Ceyhz_xn.data[ind] = -1/(dxepsdt+hdx*pex);

				/* Coefficients updating Hzx */
				Chzxh_xn.data[ind] = (2*mu_0-dt*pmx)/(2*mu_0+dt*pmx);
				Chzxey_xn.data[ind++] = -1/(dxmudt+hdx*pmx);
			}		
		}
	}

	/*====================================================================================================================*/
	if(IsPMLxp){

		sigma_max = tsigmax/(dx*PMLCellxp);

		init_coordinates(&Ceye_xp,PMLCellxp,ny,0.0);		
		init_coordinates(&Ceyhz_xp,PMLCellxp,ny,0.0);

		init_coordinates(&Chzxh_xp, Hzx_xp.nx,Hzx_xp.ny,0.0);		
		init_coordinates(&Chzxey_xp,Hzx_xp.nx,Hzx_xp.ny,0.0);

		init_coordinates(&Chzyh_xp, Hzy_xp.nx,Hzy_xp.ny,1.0);		
		init_coordinates(&Chzyex_xp,Hzy_xp.nx,Hzy_xp.ny,dt/(dy*mu_0));	

		for(ind=i=0;i<PMLCellxp;++i){
			rho_e = (i+0.25)/PMLCellxp;
			rho_m = (i+0.75)/PMLCellxp;
			pex = sigma_max*pow(rho_e,pml_order);
			pmx = mu_d_eps*sigma_max*pow(rho_m,pml_order);

			for(j=0;j<ny;j++){
				/* Coefficients updating Ey */
				Ceye_xp.data[ind] = (2*eps_0-dt*pex)/(2*eps_0+dt*pex);
				Ceyhz_xp.data[ind] = -1/(dxepsdt+hdx*pex);

				/* Coefficients updating Hzx */
				Chzxh_xp.data[ind] = (2*mu_0-dt*pmx)/(2*mu_0+dt*pmx);
				Chzxey_xp.data[ind] = -1/(dxmudt+hdx*pmx);
				ind++;
			}
		}
	}
	/*=====================================================================================================================*/
	if(IsPMLyn){

		sigma_max = tsigmax/(dy*PMLCellyn);


		init_coordinates(&Cexe_yn,nx,PMLCellyn,0.0);		
		init_coordinates(&Cexhz_yn,nx,PMLCellyn,0.0);

		init_coordinates(&Chzxh_yn, Hzx_yn.nx,Hzx_yn.ny,1.0);
		init_coordinates(&Chzxey_yn,Hzx_yn.nx,Hzx_yn.ny,-dt/(dx*mu_0));

		init_coordinates(&Chzyh_yn, Hzy_yn.nx,Hzy_yn.ny,0.0);		
		init_coordinates(&Chzyex_yn,Hzy_yn.nx,Hzy_yn.ny,0.0);


		for(j=0;j<PMLCellyn;j++){
			rho_e = 1-(j-0.75)/PMLCellyn;
			rho_m = 1-(j-0.25)/PMLCellyn;
			pex = sigma_max*pow(rho_e,pml_order);
			pmx = mu_d_eps*sigma_max*pow(rho_m,pml_order);
			ind = j;
			for(i=0;i<nx;++i){
				/* Coefficients updating Ey */
				Cexe_yn.data[ind] = (2*eps_0-dt*pex)/(2*eps_0+dt*pex);
				Cexhz_yn.data[ind] = 1/(dyepsdt+hdy*pex);

				/* Coefficients updating Hzy */
				Chzyh_yn.data[ind] = (2*mu_0-dt*pmx)/(2*mu_0+dt*pmx);
				Chzyex_yn.data[ind] = 1/(dymudt+hdy*pmx);
				ind += PMLCellyn;
			}
		}
	}
	/*=====================================================================================================================*/
	if(IsPMLyp){

		sigma_max = tsigmax/(dy*PMLCellyp);

		init_coordinates(&Cexe_yp,nx,PMLCellyp,0.0);
		init_coordinates(&Cexhz_yp,nx,PMLCellyp,0.0);

		init_coordinates(&Chzxh_yp, Hzx_yp.nx,Hzx_yp.ny,1.0);
		init_coordinates(&Chzxey_yp,Hzx_yp.nx,Hzx_yp.ny,-dt/(dx*mu_0));

		init_coordinates(&Chzyh_yp,Hzy_yp.nx,Hzy_yp.ny,0.0);
		init_coordinates(&Chzyex_yp,Hzy_yp.nx,Hzy_yp.ny,0.0);

		for(j=0;j<PMLCellyp;j++){
			rho_e = (j+0.25)/PMLCellyp;
			rho_m = (j+0.75)/PMLCellyp;
			pex = sigma_max*pow(rho_e,pml_order);
			pmx = mu_d_eps*sigma_max*pow(rho_m,pml_order);
			ind = j;
			for(i=0;i<nx;++i){
				/* Coefficients updating Ey */
				Cexe_yp.data[ind] = (2*eps_0-dt*pex)/(2*eps_0+dt*pex);
				Cexhz_yp.data[ind] = 1/(dyepsdt+hdy*pex);

				/* Coefficients updating Hzy */
				Chzyh_yp.data[ind] = (2*mu_0-dt*pmx)/(2*mu_0+dt*pmx);
				Chzyex_yp.data[ind] = 1/(dymudt+hdy*pmx);
				ind += PMLCellyp;
			}
		}
	}

	printf("End of initializing pml boundary conditions 2d TMx...\n");
}


void InitPMLBound(){

	/* determine the boundaries of the non-pml region */
	if(IsAnySidePML){
		if(IsTMx)
			InitPMLBoundTMx();
		if(IsTEx)
			InitPMLBoundTEx();
	}
}

#endif
