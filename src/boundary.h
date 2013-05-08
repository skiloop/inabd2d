#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include"defbnd.h"
#include"initbnd.h"
#include"updbnd.h"

void FreePMLSpace(){
	if(IsTMx){
	
		printf("Freeing TMx PML space...\n");
		
		free(Hzx_xn.data);
		free(Hzx_xp.data);

		free(Hzx_yn.data);
		free(Hzx_yp.data);

		free(Hzy_xn.data);
		free(Hzy_xp.data);

		free(Hzy_yn.data);
		free(Hzy_yp.data);

	/*=====================================================================================================================*/
		if(IsPMLxn){
			free(Ceye_xn.data);		
			free(Ceyhz_xn.data);

			free(Chzxh_xn.data);	
			free(Chzxey_xn.data);

			free(Chzyh_xn.data);
			free(Chzyex_xn.data);
		}

	/*====================================================================================================================*/
		if(IsPMLxp){
			free(Ceye_xp.data);		
			free(Ceyhz_xp.data);

			free(Chzxh_xp.data);		
			free(Chzxey_xp.data);

			free(Chzyh_xp.data);
			free(Chzyex_xp.data);
		}
	/*=====================================================================================================================*/
		if(IsPMLyn){
			
			free(Cexe_yn.data);		
			free(Cexhz_yn.data);

			free(Chzxh_yn.data);
			free(Chzxey_yn.data);

			free(Chzyh_yn.data);
			free(Chzyex_yn.data);

			
		}
	/*=====================================================================================================================*/
		if(IsPMLyp){

			
			free(Cexe_yp.data);
			free(Cexhz_yp.data);

			free(Chzxh_yp.data);
			free(Chzxey_yp.data);

			free(Chzyh_yp.data);
			free(Chzyex_yp.data);
			
		}
		
		printf("End of freeing TMx  PML space...\n");
	}

    if(IsTEx){
	

		printf("Freeing TEx PML space...\n");
		
		free(Ezx_xn.data);
		free(Ezx_xp.data);
		free(Ezx_yn.data);
		free(Ezx_yp.data);
		free(Ezy_xn.data);
		free(Ezy_xp.data);
		free(Ezy_yn.data);
		free(Ezy_yp.data);

		if(IsPMLxn){
			free(Chyh_xn.data);
			free(Chyez_xn.data);

			free(Cezxe_xn.data);
			free(Cezxhy_xn.data);

			free(Cezye_xn.data);
			free(Cezyhx_xn.data);
		}

		if(IsPMLyn){

			
			free(Chxh_yn.data);
			free(Chxez_yn.data);

			
			free(Cezye_yn.data);
			free(Cezyhx_yn.data);

			
			free(Cezxe_yn.data);
			free(Cezxhy_yn.data);
			
		}

		if(IsPMLxp){
			free(Chyh_xp.data);
			free(Chyez_xp.data);

			free(Cezxe_xp.data);
			free(Cezxhy_xp.data);		

			
			free(Cezye_xp.data);
			free(Cezyhx_xp.data);
		}
		if(IsPMLyp){
			free(Chxh_yp.data);
			free(Chxez_yp.data);
			
			free(Cezye_yp.data);
			free(Cezyhx_yp.data);
			
			free(Cezxe_yp.data);
			free(Cezxhy_yp.data);
			
		}
		printf("End of freeing TEx PML Space...\n");
	}
}

#endif 
