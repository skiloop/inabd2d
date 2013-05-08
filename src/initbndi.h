#ifndef _INIT_BNDI_H_
#define _INIT_BNDI_H_


void InitPMLBoundTExi(){
	int i,j,ind;
	
	MDOUBLE	*rho_e;
	MDOUBLE	*rho_m;
	
	COORDINATES	sigma_pey_yn_i;
	COORDINATES	sigma_pmy_yn_i;

	COORDINATES	sigma_pex_xp_i;
	COORDINATES	sigma_pmx_xp_i;

	COORDINATES	sigma_pex_xn_i;
	COORDINATES	sigma_pmx_xn_i;

	COORDINATES	sigma_pey_yp_i;
	COORDINATES	sigma_pmy_yp_i;

	MDOUBLE	sigma_max;

	printf("Initializing pml boundary conditions 2d TEx...\n");
	//printf("nxp1=%d\nnym1=%d\nny=%d\n",nxp1,nym1,ny);
	init_coordinates(&Ezx_xn_i,PMLCellxn,nym1,0.0);
	init_coordinates(&Ezx_xp_i,PMLCellxp,nym1,0.0);
	init_coordinates(&Ezx_yn_i,nxm1-PMLCellxn-PMLCellxp,PMLCellyn,0.0);
	init_coordinates(&Ezx_yp_i,nxm1-PMLCellxn-PMLCellxp,PMLCellyp,0.0);	
	init_coordinates(&Ezy_xn_i,PMLCellxn,nym1-PMLCellyn-PMLCellyp,0.0);
	init_coordinates(&Ezy_xp_i,PMLCellxp,nym1-PMLCellyn-PMLCellyp,0.0);
	init_coordinates(&Ezy_yn_i,nxm1,PMLCellyn,0.0);
	init_coordinates(&Ezy_yp_i,nxm1,PMLCellyp,0.0);	

	if(IsPMLxn){

		rho_e = (MDOUBLE*)malloc(PMLCellxn*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellxn*sizeof(MDOUBLE));

		init_coordinates(&sigma_pex_xn_i,PMLCellxn,nym1,0.0);
		init_coordinates(&sigma_pmx_xn_i,PMLCellxn,nym1,0.0);		
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxn);

		for(i=PMLCellxn;i>0;--i){
			rho_e[PMLCellxn-i] = (i-0.75)/PMLCellxn;
			rho_m[PMLCellxn-i] = (i-0.25)/PMLCellxn;
		}

		for(i=0;i<PMLCellxn;i++)
			for(j=0;j<nym1;++j){
				sigma_pex_xn_i.data[i*nym1+j] = sigma_max*pow(rho_e[i],pml_order);
				sigma_pmx_xn_i.data[i*nym1+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
			}


		init_coordinates(&Chyh_xn_i,PMLCellxn,nym1,0.0);
		init_coordinates(&Chyez_xn_i,PMLCellxn,nym1,0.0);

		init_coordinates(&Cezxe_xn_i,Ezx_xn_i.nx,Ezx_xn_i.ny,0.0);
		init_coordinates(&Cezxhy_xn_i,Ezx_xn_i.nx,Ezx_xn_i.ny,0.0);

		init_coordinates(&Cezye_xn_i ,Ezy_xn_i.nx,Ezy_xn_i.ny,1.0);
		init_coordinates(&Cezyhx_xn_i,Ezy_xn_i.nx,Ezy_xn_i.ny,-dt/(dy*eps_0));
		
		ind = 0;
		for(i=0;i<PMLCellxn;++i)
			for(j=0;j<nym1;j++){
				/* Coefficients updating Hx */
				Chyh_xn_i.data[ind] = (2*mu_0-dt*sigma_pmx_xn_i.data[ind])/(2*mu_0+dt*sigma_pmx_xn_i.data[ind]);
				Chyez_xn_i.data[ind] = (2*dt/dx)/(2*mu_0+dt*sigma_pmx_xn_i.data[ind]);
				
				/* Coefficients updating Ezx */
				Cezxe_xn_i.data[ind] = (2*eps_0-dt*sigma_pex_xn_i.data[ind])/(2*eps_0+dt*sigma_pex_xn_i.data[ind]);
				Cezxhy_xn_i.data[ind] = (2*dt/dx)/(2*eps_0+dt*sigma_pex_xn_i.data[ind]);

				ind++;
			}
		
		free(sigma_pmx_xn_i.data);
		free(sigma_pex_xn_i.data);

		free(rho_e);
		free(rho_m);
	}

	if(IsPMLyn){

		rho_e = (MDOUBLE*)malloc(PMLCellyn*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellyn*sizeof(MDOUBLE));

		init_coordinates(&sigma_pey_yn_i,nxm1,PMLCellyn,0.0);
		init_coordinates(&sigma_pmy_yn_i,nxm1,PMLCellyn,0.0);		
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyn);

		for(i=PMLCellyn;i>0;--i){
			rho_e[PMLCellyn-i] = (i-0.75)/PMLCellyn;
			rho_m[PMLCellyn-i] = (i-0.25)/PMLCellyn;
		}

		for(i=0;i<PMLCellyn;i++)
			for(j=0;j<nxm1;++j){
				sigma_pey_yn_i.data[j*PMLCellyn+i] = sigma_max*pow(rho_e[i],pml_order);
				sigma_pmy_yn_i.data[j*PMLCellyn+i] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
			}

		init_coordinates(&Chxh_yn_i,nxm1,PMLCellyn,0.0);
		init_coordinates(&Chxez_yn_i,nxm1,PMLCellyn,0.0);

		
		init_coordinates(&Cezye_yn_i,Ezy_yn_i.nx,Ezy_yn_i.ny,0.0);
		init_coordinates(&Cezyhx_yn_i,Ezy_yn_i.nx,Ezy_yn_i.ny,0.0);
		
		ind = 0;
		for(i=0;i<nxm1;++i)
			for(j=0;j<PMLCellyn;j++){
				/* Coefficients updating Hx */
				Chxh_yn_i.data[ind] = (2*mu_0-dt*sigma_pmy_yn_i.data[ind])/(2*mu_0+dt*sigma_pmy_yn_i.data[ind]);
				Chxez_yn_i.data[ind] = -(2*dt/dy)/(2*mu_0+dt*sigma_pmy_yn_i.data[ind]);
			
				/* Coefficients updating Ezy */
				Cezye_yn_i.data[ind] = (2*eps_0-dt*sigma_pey_yn_i.data[ind])/(2*eps_0+dt*sigma_pey_yn_i.data[ind]);
				Cezyhx_yn_i.data[ind] = -(2*dt/dy)/(2*eps_0+dt*sigma_pey_yn_i.data[ind]);
				ind++;
			}

		
		init_coordinates(&Cezxe_yn_i,Ezx_yn_i.nx,Ezx_yn_i.ny,1.0);
		init_coordinates(&Cezxhy_yn_i,Ezx_yn_i.nx,Ezx_yn_i.ny,dt/(dx*eps_0));
		
		free(sigma_pmy_yn_i.data);
		free(sigma_pey_yn_i.data);

		free(rho_e);
		free(rho_m);
	}

	if(IsPMLxp){

		rho_e = (MDOUBLE*)malloc(PMLCellxp*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellxp*sizeof(MDOUBLE));

		init_coordinates(&sigma_pex_xp_i,PMLCellxp,nym1,0.0);
		init_coordinates(&sigma_pmx_xp_i,PMLCellxp,nym1,0.0);
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxp);

		for(i=0;i<PMLCellyp;++i){
			rho_e[i] = (i+1-0.75)/PMLCellxp;
			rho_m[i] = (i+1-0.25)/PMLCellxp;
		}

		for(i=0;i<nym1;i++)
			for(j=0;j<PMLCellxp;++j){
				sigma_pex_xp_i.data[j*nym1+i] = sigma_max*pow(rho_e[j],pml_order);
				sigma_pmx_xp_i.data[j*nym1+i] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
			}
		
		init_coordinates(&Chyh_xp_i,PMLCellxp,nym1,0.0);
		init_coordinates(&Chyez_xp_i,PMLCellxp,nym1,0.0);

		init_coordinates(&Cezxe_xp_i,Ezx_xp_i.nx,Ezx_xp_i.ny,0.0);
		init_coordinates(&Cezxhy_xp_i,Ezx_xp_i.nx,Ezx_xp_i.ny,0.0);		

		ind = 0;
		for(i=0;i<PMLCellxp;++i)
			for(j=0;j<nym1;j++){
				/* Coefficients updating Hy */
				Chyh_xp_i.data[ind] = (2*mu_0-dt*sigma_pmx_xp_i.data[ind])/(2*mu_0+dt*sigma_pmx_xp_i.data[ind]);
				Chyez_xp_i.data[ind] = (2*dt/dx)/(2*mu_0+dt*sigma_pmx_xp_i.data[ind]);
				
				/* Coefficients updating Ezx */
				Cezxe_xp_i.data[ind] = (2*eps_0-dt*sigma_pex_xp_i.data[ind])/(2*eps_0+dt*sigma_pex_xp_i.data[ind]);
				Cezxhy_xp_i.data[ind] = (2*dt/dx)/(2*eps_0+dt*sigma_pex_xp_i.data[ind]);

				ind++;
			}

		
		init_coordinates(&Cezye_xp_i,Ezy_xp_i.nx,Ezy_xp_i.ny,1.0);
		init_coordinates(&Cezyhx_xp_i,Ezy_xp_i.nx,Ezy_xp_i.ny,-dt/(dy*eps_0));
		
		free(sigma_pmx_xp_i.data);
		free(sigma_pex_xp_i.data);

		free(rho_e);
		free(rho_m);
	}
	if(IsPMLyp){

		rho_e = (MDOUBLE*)malloc(PMLCellyp*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellyp*sizeof(MDOUBLE));

		init_coordinates(&sigma_pey_yp_i,nxm1,PMLCellyp,0.0);
		init_coordinates(&sigma_pmy_yp_i,nxm1,PMLCellyp,0.0);
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyp);

		for(i=0;i<PMLCellyp;++i){
			rho_e[i] = (i+1-0.75)/PMLCellyp;
			rho_m[i] = (i+1-0.25)/PMLCellyp;
		}

		for(i=0;i<nxm1;i++)
			for(j=0;j<PMLCellyp;++j){
				sigma_pey_yp_i.data[i*PMLCellyp+j] = sigma_max*pow(rho_e[j],pml_order);
				sigma_pmy_yp_i.data[i*PMLCellyp+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
			}
		
		init_coordinates(&Chxh_yp_i,nxm1,PMLCellyp,0.0);
		init_coordinates(&Chxez_yp_i,nxm1,PMLCellyp,0.0);
		
		init_coordinates(&Cezye_yp_i, Ezy_yp_i.nx,Ezy_yp_i.ny,0.0);
		init_coordinates(&Cezyhx_yp_i,Ezy_yp_i.nx,Ezy_yp_i.ny,0.0);

		ind = 0;
		for(i=0;i<nxm1;++i)
			for(j=0;j<PMLCellyp;j++){
				/* Coefficients updating Hy */
				Chxh_yp_i.data[ind] = (2*mu_0-dt*sigma_pmy_yp_i.data[ind])/(2*mu_0+dt*sigma_pmy_yp_i.data[ind]);
				Chxez_yp_i.data[ind] = -(2*dt/dx)/(2*mu_0+dt*sigma_pmy_yp_i.data[ind]);
								
				/* Coefficients updating Ezy */
				Cezye_yp_i.data[ind] = (2*eps_0-dt*sigma_pey_yp_i.data[ind])/(2*eps_0+dt*sigma_pey_yp_i.data[ind]);
				Cezyhx_yp_i.data[ind] = -(2*dt/dy)/(2*eps_0+dt*sigma_pey_yp_i.data[ind]);
				ind++;
			}
		
		init_coordinates(&Cezxe_yp_i, Ezx_yp_i.nx,Ezx_yp_i.ny,1.0);
		init_coordinates(&Cezxhy_yp_i,Ezx_yp_i.nx,Ezx_yp_i.ny,dt/(dx*eps_0));
		
		free(sigma_pmy_yp_i.data);
		free(sigma_pey_yp_i.data);

		free(rho_e);
		free(rho_m);
	}
	

	printf("End of initializing pml boundary conditions 2d TEx...\n");
}
void InitPMLBoundTMxi(){

	int i,j,ind;

	MDOUBLE	*rho_e;
	MDOUBLE	*rho_m;

	COORDINATES	sigma_pex_xn_i;
	COORDINATES	sigma_pmx_xn_i;

	COORDINATES	sigma_pey_yn_i;
	COORDINATES	sigma_pmy_yn_i;

	COORDINATES	sigma_pex_xp_i;
	COORDINATES	sigma_pmx_xp_i;

	COORDINATES	sigma_pey_yp_i;
	COORDINATES	sigma_pmy_yp_i;
	

	MDOUBLE	sigma_max;

	printf("Initializing pml boundary conditions 2d TMx...\n");
	
	init_coordinates(&Hzx_xn_i,PMLCellxn,ny,0.0);
	init_coordinates(&Hzx_xp_i,PMLCellxp,ny,0.0);

	init_coordinates(&Hzx_yn_i,nx-PMLCellxn-PMLCellxp,PMLCellyn,0.0);
	init_coordinates(&Hzx_yp_i,nx-PMLCellxp-PMLCellxn,PMLCellyp,0.0);

	init_coordinates(&Hzy_xn_i,PMLCellxn,ny-PMLCellyn-PMLCellyp,0.0);
	init_coordinates(&Hzy_xp_i,PMLCellxp,ny-PMLCellyp-PMLCellyn,0.0);

	init_coordinates(&Hzy_yn_i,nx,PMLCellyn,0.0);
	init_coordinates(&Hzy_yp_i,nx,PMLCellyp,0.0);

/*=====================================================================================================================*/
	if(IsPMLxn){
		rho_e = (MDOUBLE*)malloc(PMLCellxn*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellxn*sizeof(MDOUBLE));		

		init_coordinates(&sigma_pex_xn_i,PMLCellxn,ny,0.0);		
		init_coordinates(&sigma_pmx_xn_i,PMLCellxn,ny,0.0);
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxn);

		for(i=PMLCellxn;i>0;--i){
			rho_e[PMLCellxn-i] = (i-0.75)/PMLCellxn;
			rho_m[PMLCellxn-i] = (i-0.25)/PMLCellxn;
		}

		for(i=0;i<PMLCellxn;i++)
			for(j=0;j<ny;++j){
				sigma_pex_xn_i.data[i*ny+j] = sigma_max*pow(rho_e[i],pml_order);
				sigma_pmx_xn_i.data[i*ny+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
			}

		init_coordinates(&Ceye_xn_i,PMLCellxn,ny,0.0);		
		init_coordinates(&Ceyhz_xn_i,PMLCellxn,ny,0.0);

		init_coordinates(&Chzxh_xn_i, Hzx_xn_i.nx,Hzx_xn_i.ny,0.0);		
		init_coordinates(&Chzxey_xn_i,Hzx_xn_i.nx,Hzx_xn_i.ny,0.0);

		init_coordinates(&Chzyh_xn_i, Hzy_xn_i.nx,Hzy_xn_i.ny,1.0);	
		init_coordinates(&Chzyex_xn_i,Hzy_xn_i.nx,Hzy_xn_i.ny,dt/(dy*mu_0));

		ind = 0;
		for(i=0;i<PMLCellxn;++i)
			for(j=0;j<ny;j++){
				/* Coefficients updating Ey */
				Ceye_xn_i.data[ind] = (2*eps_0-dt*sigma_pex_xn_i.data[ind])/(2*eps_0+dt*sigma_pex_xn_i.data[ind]);
				Ceyhz_xn_i.data[ind] = -(2*dt/dx)/(2*eps_0+dt*sigma_pex_xn_i.data[ind]);
			
				/* Coefficients updating Hzx */
				Chzxh_xn_i.data[ind] = (2*mu_0-dt*sigma_pmx_xn_i.data[ind])/(2*mu_0+dt*sigma_pmx_xn_i.data[ind]);
				Chzxey_xn_i.data[ind] = -(2*dt/dx)/(2*mu_0+dt*sigma_pmx_xn_i.data[ind]);

				ind++;
			}
		
		free(sigma_pmx_xn_i.data);
		free(sigma_pex_xn_i.data);

		free(rho_e);
		free(rho_m);
	}

/*====================================================================================================================*/
	if(IsPMLxp){

		rho_e = (MDOUBLE*)malloc(PMLCellxp*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellxp*sizeof(MDOUBLE));

		init_coordinates(&sigma_pex_xp_i,PMLCellxp,ny,0.0);		
		init_coordinates(&sigma_pmx_xp_i,PMLCellxp,ny,0.0);
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dx*PMLCellxp);

		for(i=0;i<PMLCellxp;++i){
			rho_e[i] = (i+1-0.75)/PMLCellxp;
			rho_m[i] = (i+1-0.25)/PMLCellxp;
		}

		for(i=0;i<PMLCellxp;i++)
			for(j=0;j<ny;++j){
				sigma_pex_xp_i.data[i*ny+j] = sigma_max*pow(rho_e[i],pml_order);
				sigma_pmx_xp_i.data[i*ny+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[i],pml_order);
			}

		init_coordinates(&Ceye_xp_i,PMLCellxp,ny,0.0);		
		init_coordinates(&Ceyhz_xp_i,PMLCellxp,ny,0.0);

		init_coordinates(&Chzxh_xp_i, Hzx_xp_i.nx,Hzx_xp_i.ny,0.0);		
		init_coordinates(&Chzxey_xp_i,Hzx_xp_i.nx,Hzx_xp_i.ny,0.0);

		init_coordinates(&Chzyh_xp_i, Hzy_xp_i.nx,Hzy_xp_i.ny,1.0);		
		init_coordinates(&Chzyex_xp_i,Hzy_xp_i.nx,Hzy_xp_i.ny,dt/(dy*mu_0));	

		ind = 0;
		for(i=0;i<PMLCellxp;++i)
			for(j=0;j<ny;j++){
				/* Coefficients updating Ey */
				Ceye_xp_i.data[ind] = (2*eps_0-dt*sigma_pex_xp_i.data[ind])/(2*eps_0+dt*sigma_pex_xp_i.data[ind]);
				Ceyhz_xp_i.data[ind] = -(2*dt/dx)/(2*eps_0+dt*sigma_pex_xp_i.data[ind]);
				
				/* Coefficients updating Hzx */
				Chzxh_xp_i.data[ind] = (2*mu_0-dt*sigma_pmx_xp_i.data[ind])/(2*mu_0+dt*sigma_pmx_xp_i.data[ind]);
				Chzxey_xp_i.data[ind] = -(2*dt/dx)/(2*mu_0+dt*sigma_pmx_xp_i.data[ind]);

				ind++;
			}
		
		free(sigma_pmx_xp_i.data);
		free(sigma_pex_xp_i.data);
		
		free(rho_e);
		free(rho_m);
		
	}
/*=====================================================================================================================*/
	if(IsPMLyn){
		rho_e = (MDOUBLE*)malloc(PMLCellyn*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellyn*sizeof(MDOUBLE));

		init_coordinates(&sigma_pey_yn_i,nx,PMLCellyn,0.0);		
		init_coordinates(&sigma_pmy_yn_i,nx,PMLCellyn,0.0);
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyn);

		for(i=PMLCellyn;i>0;--i){
			rho_e[PMLCellyn-i] = (i-0.75)/PMLCellyn;
			rho_m[PMLCellyn-i] = (i-0.25)/PMLCellyn;
		}

		for(i=0;i<nx;i++)
			for(j=0;j<PMLCellyn;++j){
				sigma_pey_yn_i.data[i*PMLCellyn+j] = sigma_max*pow(rho_e[j],pml_order);
				sigma_pmy_yn_i.data[i*PMLCellyn+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
			}

		init_coordinates(&Cexe_yn_i,nx,PMLCellyn,0.0);		
		init_coordinates(&Cexhz_yn_i,nx,PMLCellyn,0.0);

		init_coordinates(&Chzxh_yn_i, Hzx_yn_i.nx,Hzx_yn_i.ny,1.0);
		init_coordinates(&Chzxey_yn_i,Hzx_yn_i.nx,Hzx_yn_i.ny,-dt/(dx*mu_0));

		init_coordinates(&Chzyh_yn_i, Hzy_yn_i.nx,Hzy_yn_i.ny,0.0);		
		init_coordinates(&Chzyex_yn_i,Hzy_yn_i.nx,Hzy_yn_i.ny,0.0);

		ind = 0;
		for(i=0;i<nx;++i)
			for(j=0;j<PMLCellyn;j++){
				/* Coefficients updating Ey */
				Cexe_yn_i.data[ind] = (2*eps_0-dt*sigma_pey_yn_i.data[ind])/(2*eps_0+dt*sigma_pey_yn_i.data[ind]);
				Cexhz_yn_i.data[ind] = (2*dt/dy)/(2*eps_0+dt*sigma_pey_yn_i.data[ind]);
				
				/* Coefficients updating Hzy */
				Chzyh_yn_i.data[ind] = (2*mu_0-dt*sigma_pmy_yn_i.data[ind])/(2*mu_0+dt*sigma_pmy_yn_i.data[ind]);
				Chzyex_yn_i.data[ind] = (2*dt/dy)/(2*mu_0+dt*sigma_pmy_yn_i.data[ind]);
				ind++;
			}
		
		free(sigma_pmy_yn_i.data);
		free(sigma_pey_yn_i.data);

		free(rho_e);
		free(rho_m);
	}
/*=====================================================================================================================*/
	if(IsPMLyp){

		rho_e = (MDOUBLE*)malloc(PMLCellyp*sizeof(MDOUBLE));
		rho_m = (MDOUBLE*)malloc(PMLCellyp*sizeof(MDOUBLE));

		init_coordinates(&sigma_pey_yp_i,nx,PMLCellyp,0.0);
		init_coordinates(&sigma_pmy_yp_i,nx,PMLCellyp,0.0);
		
		sigma_max = -(pml_order+1)*eps_0*c*log(pml_R_0)/(2*dy*PMLCellyp);

		for(i=0;i<PMLCellyp;++i){
			rho_e[i] = (i+1-0.75)/PMLCellyp;
			rho_m[i] = (i+1-0.25)/PMLCellyp;
		}

		for(i=0;i<nx;i++)
			for(j=0;j<PMLCellyp;++j){
				sigma_pey_yp_i.data[i*PMLCellyp+j] = sigma_max*pow(rho_e[j],pml_order);
				sigma_pmy_yp_i.data[i*PMLCellyp+j] = (mu_0/eps_0)*sigma_max*pow(rho_m[j],pml_order);
			}

		init_coordinates(&Cexe_yp_i,nx,PMLCellyp,0.0);
		init_coordinates(&Cexhz_yp_i,nx,PMLCellyp,0.0);

		init_coordinates(&Chzxh_yp_i, Hzx_yp_i.nx,Hzx_yp_i.ny,1.0);
		init_coordinates(&Chzxey_yp_i,Hzx_yp_i.nx,Hzx_yp_i.ny,-dt/(dx*mu_0));

		init_coordinates(&Chzyh_yp_i,Hzy_yp_i.nx,Hzy_yp_i.ny,0.0);
		init_coordinates(&Chzyex_yp_i,Hzy_yp_i.nx,Hzy_yp_i.ny,0.0);
		ind = 0;
		for(i=0;i<nx;++i)
			for(j=0;j<PMLCellyp;j++){
				/* Coefficients updating Ey */
				Cexe_yp_i.data[ind] = (2*eps_0-dt*sigma_pey_yp_i.data[ind])/(2*eps_0+dt*sigma_pey_yp_i.data[ind]);
				Cexhz_yp_i.data[ind] = (2*dt/dy)/(2*eps_0+dt*sigma_pey_yp_i.data[ind]);
								
				/* Coefficients updating Hzy */
				Chzyh_yp_i.data[ind] = (2*mu_0-dt*sigma_pmy_yp_i.data[ind])/(2*mu_0+dt*sigma_pmy_yp_i.data[ind]);
				Chzyex_yp_i.data[ind] = (2*dt/dy)/(2*mu_0+dt*sigma_pmy_yp_i.data[ind]);
				ind++;
			}
	
	free(sigma_pmy_yp_i.data);
	free(sigma_pey_yp_i.data);

	free(rho_e);
	free(rho_m);
	}
	
	printf("End of initializing pml boundary conditions 2d TMx...\n");
}


void InitPMLBoundi(){

	/* determine the boundaries of the non-pml region */
	if(IsAnySidePML){
		if(IsTMx)
			InitPMLBoundTMxi();
		if(IsTEx)
			InitPMLBoundTExi();
	}
}




#endif
