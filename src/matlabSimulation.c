
#include "common.h"
#include "commonData.h"

#ifdef MATLAB_SIMULATION

//#include"matrix.h"
#include"engine.h"
#include"mex.h"
#ifdef printf
#undef printf
#endif

Engine *ep = NULL;

#include "dataType.h"

void PreSimuAdd1(Engine *ep, const MyStruct stru) {
    mxArray *MyArray;
    CheckStruct(stru);
    MyArray = mxCreateDoubleMatrix(stru.ny, stru.nx, mxREAL);
    memcpy(mxGetPr(MyArray), stru.data, stru.nx * stru.ny * sizeof (MyDataF));
    engPutVariable(ep, "data1", MyArray);
    engEvalString(ep, "figure,h1=imagesc(data1);set(gca,'YDir','normal');"); //surf(data1);");//
    engEvalString(ep, "axes_h1 = gca;");
    engEvalString(ep, "ind1=1;xlabel(axes_h1,ind1);");
    engEvalString(ep, "colorbar");
    engEvalString(ep, "title('Display Data1');");
    engEvalString(ep, "set(gcf,'Renderer','zbuffer');");
    engEvalString(ep, "hold on");
    mxDestroyArray(MyArray);

}

void PreSimuAdd0(Engine *ep, const MyStruct stru) {

    mxArray *MyArray;
    CheckStruct(stru);
    MyArray = mxCreateDoubleMatrix(stru.ny, stru.nx, mxREAL);
    memcpy(mxGetPr(MyArray), stru.data, stru.nx * stru.ny * sizeof (MyDataF));
    engPutVariable(ep, "data0", MyArray);
    engEvalString(ep, "figure,h0=imagesc(data0);set(gca,'YDir','normal');");
    engEvalString(ep, "axes_h0 = gca;");
    engEvalString(ep, "ind0=1;xlabel(axes_h0,ind0);");
    engEvalString(ep, "colorbar");
    engEvalString(ep, "title('Display Data0');");
    engEvalString(ep, "set(gcf,'Renderer','zbuffer');");
    engEvalString(ep, "hold on");
    mxDestroyArray(MyArray);

}

void PreSimuAdd2(Engine *ep, const MyStruct stru) {
    mxArray *MyArray;
    CheckStruct(stru);
    MyArray = mxCreateDoubleMatrix(stru.ny, stru.nx, mxREAL);
    memcpy(mxGetPr(MyArray), stru.data, stru.nx * stru.ny * sizeof (MyDataF));
    engPutVariable(ep, "data2", MyArray);
    engEvalString(ep, "figure,h2=imagesc(data2);set(gca,'YDir','normal');");
    engEvalString(ep, "axes_h2 = gca;");
    engEvalString(ep, "ind2=1;xlabel(axes_h2,ind2);");
    engEvalString(ep, "colorbar");
    engEvalString(ep, "title('Display Data2');");
    engEvalString(ep, "set(gcf,'Renderer','zbuffer');");
    engEvalString(ep, "hold on");
    mxDestroyArray(MyArray);

}

void Simulate0(Engine *ep, const MyStruct stru) {
    mxArray *MyArray;
    CheckStruct(stru);
    MyArray = mxCreateDoubleMatrix(stru.ny, stru.nx, mxREAL);
    memcpy(mxGetPr(MyArray), stru.data, stru.nx * stru.ny * sizeof (MyDataF));
    engPutVariable(ep, "data0", MyArray);
    engEvalString(ep, "ind0 = ind0 + 1;");
    engEvalString(ep, "set(h0,'cdata',data0);");
    engEvalString(ep, "xlabel(axes_h0,ind0);");
    engEvalString(ep, "drawnow");
    mxDestroyArray(MyArray);
}

void Simulate1(Engine *ep, const MyStruct stru) {
    mxArray *MyArray;
    CheckStruct(stru);
    MyArray = mxCreateDoubleMatrix(stru.ny, stru.nx, mxREAL);
    memcpy(mxGetPr(MyArray), stru.data, stru.nx * stru.ny * sizeof (MyDataF));
    engPutVariable(ep, "data1", MyArray);
    engEvalString(ep, "ind1 = ind1 + 1;");
    engEvalString(ep, "set(h1,'cdata',data1);"); //
    //engEvalString(ep,"set(h1,'zdata',data1);");//
    engEvalString(ep, "xlabel(axes_h1,ind1);");
    engEvalString(ep, "drawnow");
    mxDestroyArray(MyArray);
}

void Simulate2(Engine *ep, const MyStruct stru) {
    mxArray *MyArray;
    CheckStruct(stru);
    MyArray = mxCreateDoubleMatrix(stru.ny, stru.nx, mxREAL);
    memcpy(mxGetPr(MyArray), stru.data, stru.nx * stru.ny * sizeof (MyDataF));
    engPutVariable(ep, "data2", MyArray);
    engEvalString(ep, "ind2 = ind2 + 1;");
    engEvalString(ep, "set(h2,'cdata',data2);");
    engEvalString(ep, "xlabel(axes_h2,ind2);");
    engEvalString(ep, "drawnow");
    mxDestroyArray(MyArray);
}

void PreSimuAddNE(Engine *ep, const MyStruct stru) {
    mxArray *MyArray;
    CheckStruct(stru);
    MyArray = mxCreateDoubleMatrix(stru.ny, stru.nx, mxREAL);
    memcpy(mxGetPr(MyArray), stru.data, stru.nx * stru.ny * sizeof (MyDataF));
    engPutVariable(ep, "NE", MyArray);
    engEvalString(ep, "figure,NE_h=imagesc(NE);set(gca,'YDir','normal');");
    engEvalString(ep, "NE_axes_h = gca;");
    engEvalString(ep, "NE_ind=1;xlabel(NE_axes_h,NE_ind);");
    engEvalString(ep, "colorbar");
    engEvalString(ep, "title('Display NE');");
    engEvalString(ep, "set(gcf,'Renderer','zbuffer');");
    engEvalString(ep, "hold on");
    mxDestroyArray(MyArray);

}

void SimulateNe(Engine *ep, const MyStruct stru) {
    mxArray *MyArray;
    CheckStruct(stru);
    MyArray = mxCreateDoubleMatrix(stru.ny, stru.nx, mxREAL);
    memcpy(mxGetPr(MyArray), stru.data, stru.nx * stru.ny * sizeof (MyDataF));
    engPutVariable(ep, "NE", MyArray);
    engEvalString(ep, "NE_ind = NE_ind + 1;");
    engEvalString(ep, "set(NE_h,'cdata',NE);");
    engEvalString(ep, "xlabel(NE_axes_h,NE_ind);");
    engEvalString(ep, "drawnow");
    mxDestroyArray(MyArray);
}

#endif //if define MATLAB_SIMULTION

void DispEMFields(const int timestep) {

    if (IsMatlabSim) {
#ifdef MATLAB_SIMULATION
        //Matlab Simulation
        if (timestep % LEAPSTEP_OF_DISPLAY == 0) {
            if (IsTMx) {

#ifdef FIELD_TO_DISPLAY0
                //display 0
                Simulate0(ep, FIELD_TO_DISPLAY0);
#endif
                //display 1
#ifdef FIELD_TO_DISPLAY1
                Simulate1(ep, FIELD_TO_DISPLAY1);

#endif
                //display 2
#ifdef FIELD_TO_DISPLAY2
                Simulate2(ep, FIELD_TO_DISPLAY2);
#endif
            }
            if (IsTEx) {
#ifdef FIELD_TO_DISPLAY3
                //display 0
                Simulate0(ep, FIELD_TO_DISPLAY3);
#endif
                //display 1
#ifdef FIELD_TO_DISPLAY4
                Simulate1(ep, FIELD_TO_DISPLAY4);

#endif
                //display 2
#ifdef FIELD_TO_DISPLAY5
                Simulate2(ep, FIELD_TO_DISPLAY5);
#endif
            }
        }
#endif
    }

}

void EndMatlabSim() {
    if (IsMatlabSim) {
#ifdef MATLAB_SIMULATION
        //Close Matlab Engine
        engEvalString(ep, "clear;close all;");
        engClose(ep);
#endif
    }

}

void DispNe(const int timestep) {
    if (IsMatlabSim) {
#ifdef MATLAB_SIMULATION

        if (IfWithDensity)SimulateNe(ep, ne);

#endif	
    }

}

void InitSim() {
    if (IsMatlabSim) {
#ifdef MATLAB_SIMULATION
        //Define and open matlab engine
        if ((ep = engOpen(NULL)) == NULL) {
            printf("Can't start matlab engine!\n");
            exit(1);
        }
        if (IsTMx) {
#ifdef FIELD_TO_DISPLAY0
            PreSimuAdd0(ep, FIELD_TO_DISPLAY0);
#endif
#ifdef FIELD_TO_DISPLAY1
            PreSimuAdd1(ep, FIELD_TO_DISPLAY1);
#endif
#ifdef FIELD_TO_DISPLAY2
            PreSimuAdd2(ep, FIELD_TO_DISPLAY2);
#endif
        }
        if (IsTEx) {
#ifdef FIELD_TO_DISPLAY3
            PreSimuAdd0(ep, FIELD_TO_DISPLAY3);
#endif
#ifdef FIELD_TO_DISPLAY4
            PreSimuAdd1(ep, FIELD_TO_DISPLAY4);
#endif
#ifdef FIELD_TO_DISPLAY5
            PreSimuAdd2(ep, FIELD_TO_DISPLAY5);
#endif
        }
#ifdef DISPLAY_NE
        if (IfWithDensity)PreSimuAddNE(ep, ne);
#endif

#endif
    }
}
