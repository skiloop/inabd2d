
#include <math.h>

#include "fdtd.h"
#include "fluid.h"
#include "commonData.h"

//=========================================================================================================================

MyDataF ElectronTemperature(MyDataF Em, MyDataF p) {

    /*
	
    MyDataF Ngas = 2.44e25*p/760.0;
    MyDataF EmDivNgas = Em / Ngas *1e21;
	
    MyDataF y, y0, xc, w, A;
	
    y0 = 51.54225; 
    xc = -6456.19545;
	
    w  = 3645.62175;
    A  = -3.9399E6;
	
    y = y0 + (2*A/M_PI)*(w/(4*pow((EmDivNgas-xc), 2) + w*w));
	
    return (2.0/3.0*y);
     */
    return 2.0;

}
//=========================================================================================================================

MyDataF collision(MyDataF Em, MyDataF p) {

    /*
    MyDataF Ngas = 2.44e25*p/760.0;
    MyDataF EmDivNgas = Em / Ngas *1e21;
	
    MyDataF y, y0, xc, w, A;
	
    y0 = 3.88827E-13; 
    xc = -3042.82541;
	
    w  = 148.94905;
    A  = -1.1716E-7;
	
    y = y0 + (2*A/M_PI)*(w/(4*pow((EmDivNgas-xc), 2) + w*w));
	
    return (y*Ngas);
     */
    return 4.0e12;

}
//=======================================================================================================
/*MyDataF ionization(MyDataF Em, MyDataF p)
{
        Em = Em/100.0;
	

        if(Em/p<54) return 0.0;
        else if(Em/p<=120) return p*(5.0+0.19*(Em/p))*1.0e7*exp(-273.8*p/Em);
        else if(Em/p<3000) return p*54.08*1.0e6*sqrt(Em/p)*exp(-359*p/Em);
        else return 0.0;
}*/

//=========================================================================================================================

MyDataF ionization(MyDataF Em, MyDataF p) {
    MyDataF Ngas = 2.44e25 * p / 760.0;
    MyDataF EmDivNgas = Em / Ngas * 1e21;

    MyDataF y, y0, xc, w, A;

    if (EmDivNgas < 19.53) return 0.0;

    else if (EmDivNgas < 39.06) {
        return Ngas * (2.5935E-29 + (9.54368E-23 - 2.5935E-29) / (39.06 - 19.53)*(EmDivNgas - 19.53));
    } else if (EmDivNgas < 58.59) {
        return Ngas * (9.54368E-23 + (1.58237E-20 - 9.54368E-23) / (58.59 - 39.06)*(EmDivNgas - 39.06));
    }//if(EmDivNgas<58.59) return 0.0;

    else if (EmDivNgas < 78.13) {
        return Ngas * (1.58237E-20 + (2.49102E-19 - 1.58237E-20) / (78.13 - 58.59)*(EmDivNgas - 58.59));
    } else if (EmDivNgas < 97.66) {
        return Ngas * (2.49102E-19 + (1.44285E-18 - 2.49102E-19) / (97.66 - 78.13)*(EmDivNgas - 78.13));
    } else if (EmDivNgas < 117.2) {
        return Ngas * (1.44285E-18 + (4.98158E-18 - 1.44285E-18) / (117.2 - 97.66)*(EmDivNgas - 97.66));
    } else if (EmDivNgas < 273.4) {
        y0 = -4.61936E-17;
        xc = 293.32953;
        w = 131.54322;
        A = 8.25101E-14;
    } else if (EmDivNgas < 625) {
        y0 = -1.49081E-15;
        xc = 769.51501;
        w = 595.80256;
        A = 6.33662E-12;
    } else {
        y0 = 1.82556E-13;
        xc = -1422.48023;

        w = 11138.83498;
        A = -3.56162E-9;
    }


    y = y0 + (2 * A / M_PI)*(w / (4 * pow((EmDivNgas - xc), 2) + w * w));

    return (y * Ngas);

}


//=========================================================================================================================

MyDataF attachment(MyDataF Em, MyDataF p) {

    /*

    MyDataF Ngas = 2.44e25*p/760.0;
    MyDataF EmDivNgas = Em / Ngas *1e21;
	
    MyDataF y, y0, xc, w, A;
     */
    return 0.0;

    /*
        if(EmDivNgas<2.4) return 0.0;
        else if(EmDivNgas<9.766)
        {
                y0 = 1.51326E-37; 
        xc = 2.441;
	
            w  = 1.16603;
            A  = 3.1722E-37;
        }
   else if(EmDivNgas<19.53) return Ngas*(-1.47072E-23 + 1.50596E-24*EmDivNgas); 

   

   else if(EmDivNgas<39.06) return Ngas*( -4.47216E-20+ 2.29064E-21*EmDivNgas); 
     */
    /*
    if(EmDivNgas<39.06) return 0.0;

    else if(EmDivNgas<58.59) return Ngas*( -1.05687E-18+ 2.82032E-20*EmDivNgas); 
   
   // if(EmDivNgas<58.59) return 0.0;

    else if(EmDivNgas< 117.2) return Ngas*( -4.15977E-18+ 8.00397E-20*EmDivNgas);

    else if(EmDivNgas<351.6) 
    {
            y0 = 1.00104E-17; 
            xc =88.55566;
	   
            w  = 119.43221;
            A  = -1.09822E-15;	   
         }
    else if(EmDivNgas<2070) 
    {
            y0 = 6.93405E-18; 
            xc =283.27689;
	   
            w  = 1113.47354;
            A  = 4.80788E-15;	   
         }
   
    else 
         {
            y0 = 6.1407E-18; 
            xc = 6940.48841;	
            w  = 6125.1287;
            A  = 3.44051E-14;
         }
   
    y = y0 + (2*A/M_PI)*(w/(4*pow((EmDivNgas-xc), 2) + w*w));
    return (y*Ngas);	
     */
}



