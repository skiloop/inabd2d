#ifndef __INC_BOUND_H__
#define __INC_BOUND_H__


#include"defbndi.h"
#include"initbndi.h"
#include"updbnd.h"
void FreePMLSpacei(){
	if(IsTMx){
	
		printf("Freeing TMx PML space...\n");
		
		free(Hzx_xn_i.data);
		free(Hzx_xp_i.data);

		free(Hzx_yn_i.data);
		free(Hzx_yp_i.data);

		free(Hzy_xn_i.data);
		free(Hzy_xp_i.data);

		free(Hzy_yn_i.data);
		free(Hzy_yp_i.data);

	/*=====================================================================================================================*/
		if(IsPMLxn){
			free(Ceye_xn_i.data);		
			free(Ceyhz_xn_i.data);

			free(Chzxh_xn_i.data);	
			free(Chzxey_xn_i.data);

			free(Chzyh_xn_i.data);
			free(Chzyex_xn_i.data);
		}

	/*====================================================================================================================*/
		if(IsPMLxp){
			free(Ceye_xp_i.data);		
			free(Ceyhz_xp_i.data);

			free(Chzxh_xp_i.data);		
			free(Chzxey_xp_i.data);

			free(Chzyh_xp_i.data);
			free(Chzyex_xp_i.data);
		}
	/*=====================================================================================================================*/
		if(IsPMLyn){
			
			free(Cexe_yn_i.data);		
			free(Cexhz_yn_i.data);

			free(Chzxh_yn_i.data);
			free(Chzxey_yn_i.data);

			free(Chzyh_yn_i.data);
			free(Chzyex_yn_i.data);

			
		}
	/*=====================================================================================================================*/
		if(IsPMLyp){

			
			free(Cexe_yp_i.data);
			free(Cexhz_yp_i.data);

			free(Chzxh_yp_i.data);
			free(Chzxey_yp_i.data);

			free(Chzyh_yp_i.data);
			free(Chzyex_yp_i.data);
			
		}
		
		printf("End of freeing TMx  PML space...\n");
	}

    if(IsTEx){
	

		printf("Freeing TEx PML space...\n");
		
		free(Ezx_xn_i.data);
		free(Ezx_xp_i.data);
		free(Ezx_yn_i.data);
		free(Ezx_yp_i.data);
		free(Ezy_xn_i.data);
		free(Ezy_xp_i.data);
		free(Ezy_yn_i.data);
		free(Ezy_yp_i.data);

		if(IsPMLxn){
			free(Chyh_xn_i.data);
			free(Chyez_xn_i.data);

			free(Cezxe_xn_i.data);
			free(Cezxhy_xn_i.data);

			free(Cezye_xn_i.data);
			free(Cezyhx_xn_i.data);
		}

		if(IsPMLyn){

			
			free(Chxh_yn_i.data);
			free(Chxez_yn_i.data);

			
			free(Cezye_yn_i.data);
			free(Cezyhx_yn_i.data);

			
			free(Cezxe_yn_i.data);
			free(Cezxhy_yn_i.data);
			
		}

		if(IsPMLxp){
			free(Chyh_xp_i.data);
			free(Chyez_xp_i.data);

			free(Cezxe_xp_i.data);
			free(Cezxhy_xp_i.data);		

			
			free(Cezye_xp_i.data);
			free(Cezyhx_xp_i.data);
		}
		if(IsPMLyp){
			free(Chxh_yp_i.data);
			free(Chxez_yp_i.data);
			
			free(Cezye_yp_i.data);
			free(Cezyhx_yp_i.data);
			
			free(Cezxe_yp_i.data);
			free(Cezxhy_yp_i.data);
			
		}
		printf("End of freeing TEx PML Space...\n");
	}
}


#endif
