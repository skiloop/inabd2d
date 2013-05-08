
#ifndef _SIMULATE_H_
#define _SIMULATE_H_

#include"matlabsim.h"

#ifdef MATLAB_SIMULATION
Engine *ep=NULL;
#endif

void DispEMFields(const int timestep){

	if(IsMatlabSim){
#ifdef MATLAB_SIMULATION
		//Matlab Simulation
		if(timestep%LEAPSTEP_OF_DISPLAY==0)
		{
			if(IsTMx){
			
#ifdef FIELD_TO_DISPLAY0
				//display 0
				Simulate0(ep,FIELD_TO_DISPLAY0);
#endif
//display 1
#ifdef FIELD_TO_DISPLAY1
				Simulate1(ep,FIELD_TO_DISPLAY1);

#endif
//display 2
#ifdef FIELD_TO_DISPLAY2
				Simulate2(ep,FIELD_TO_DISPLAY2);
#endif
			}
			if(IsTEx){
#ifdef FIELD_TO_DISPLAY3
				//display 0
				Simulate0(ep,FIELD_TO_DISPLAY3);
#endif
//display 1
#ifdef FIELD_TO_DISPLAY4
				Simulate1(ep,FIELD_TO_DISPLAY4);

#endif
//display 2
#ifdef FIELD_TO_DISPLAY5
				Simulate2(ep,FIELD_TO_DISPLAY5);
#endif
			}
		}
#endif
	}

}
void EndMatlabSim(){
	if(IsMatlabSim){
#ifdef MATLAB_SIMULATION
		//Close Matlab Engine
		engEvalString(ep,"clear;close all;");
		engClose(ep);
#endif
	}

}
void DispNe(const int timestep){
	if(IsMatlabSim){	
#ifdef MATLAB_SIMULATION
	
		if(IfWithDensity)SimulateNe(ep,ne);

#endif	
	}

}
void InitSim(){
	if(IsMatlabSim){
#ifdef MATLAB_SIMULATION
	//Define and open matlab engine
	if((ep=engOpen(NULL))==NULL){
        	printf("Can't start matlab engine!\n");
        	exit(1);
    	}
	if(IsTMx){
#ifdef FIELD_TO_DISPLAY0
		PreSimuAdd0(ep,FIELD_TO_DISPLAY0);
#endif
#ifdef FIELD_TO_DISPLAY1
		PreSimuAdd1(ep,FIELD_TO_DISPLAY1);
#endif
#ifdef FIELD_TO_DISPLAY2
		PreSimuAdd2(ep,FIELD_TO_DISPLAY2);
#endif
	}
	if(IsTEx){
#ifdef FIELD_TO_DISPLAY3
		PreSimuAdd0(ep,FIELD_TO_DISPLAY3);
#endif
#ifdef FIELD_TO_DISPLAY4
		PreSimuAdd1(ep,FIELD_TO_DISPLAY4);
#endif
#ifdef FIELD_TO_DISPLAY5
		PreSimuAdd2(ep,FIELD_TO_DISPLAY5);
#endif
	}
#ifdef DISPLAY_NE
	if(IfWithDensity)PreSimuAddNE(ep,ne);
#endif
	
#endif
	}
}

#endif //ifdef _SIMULATE_H_
	
