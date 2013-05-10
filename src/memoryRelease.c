
#include "common.h"
#include "dataType.h"
#include "commonData.h"

void FreeSpace() {
    if (IsTMx) {

        //free(Ex.data);
        //free(Ey.data);
        //free(Hz.data);
        /*----------------------------------*/

        free(Ex_s.data);
        free(Ex_s_pre.data);
        //free(Ex_i.data);
        //free(Ex_i_pre.data);

        free(Vex.data);
        free(Cevx.data);
        free(Ceex.data);
        free(Cehx.data);
        /*----------------------------------*/


        /*----------------------------------*/

        free(Ey_s.data);
        free(Ey_s_pre.data);
        //free(Ey_i.data);
        //free(Ey_i_pre.data);

        free(Vey.data);
        free(Cevy.data);
        free(Ceey.data);
        free(Cehy.data);
        /*----------------------------------*/


        /*----------------------------------*/

        free(Hz_s.data);
        //free(Hz_s_pre.data);
        //free(Hz_i.data);
        //free(Hz_i_pre.data);
        /*----------------------------------*/
    }
    if (IsTEx) {
        free(Hx_i.data);
        free(Hy_i.data);
        free(Ez_i.data);

        //free(Hx.data);
        //free(Hy.data);
        //free(Ez.data);

        /*---------------< X >-------------------*/
        free(Hx_s.data);
        //free(Hx_s_pre.data);
        //free(Hx_i.data);
        //free(Hx_i_pre.data);
        /*----------------------------------*/

        /*---------------< Y >-----------------*/
        free(Hy_s.data);
        //free(Hy_s_pre.data);


        /*----------------------------------*/


        /*--------------< Z >-----------------*/
        free(Ez_s.data);
        free(Ez_s_pre.data);
        //free(Ez_i_pre.data);

        free(Vez.data);
        free(Ceez.data);
        free(Cevz.data);
        free(Cehz.data);
        /*----------------------------------*/

    }
	if(Cvxvx.data!=NULL){
		free(Cvxvx.data);
	}
	if(Cvyvy.data!=NULL){
		free(Cvyvy.data);
	}
	if(Cvzvz.data!=NULL){
		free(Cvzvz.data);
	}
	if(Cvxex.data!=NULL){
		free(Cvxex.data);
	}
	if(Cvyey.data!=NULL){
		free(Cvyey.data);
	}
	if(Cvzez.data!=NULL){
		free(Cvzez.data);
	}
    free(ne.data);
    free(Erms.data);
    free(ne_pre.data);
    free(beta.data);
	if(Nu_c.data!=NULL){
		free(Nu_c.data);
	}
}

