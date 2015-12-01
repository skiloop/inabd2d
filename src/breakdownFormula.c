

#include <math.h>

#include "common.h"
#include "dataType.h"
#include "breakdownFormula.h"

//////////////////////////////////////////////////////////////////////////
//int sign(MyDataF val)
//if val>0 return 1
//if val<0 return -1
//return 0 if val==0

int sign(MyDataF val) {
    return (val < 0 ? -1 : (val > 0 ? 1 : 0));
}
//////////////////////////////////////////////////////////////////////////
// Nikonov formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Nikonov formula

MyDataF Alpha_Nikonov(MyDataF E, MyDataF P) {
    MyDataF EDivP = fabs(E) / P;
    if (EDivP < 108.0) {
        return 3.9 * P * exp(-213.0 / EDivP);
    } else {
        return 14.5 * P * exp(-356.0 / EDivP);
    }
}
// Calculate Eta by Nikonov formula

MyDataF Eta_Nikonov(MyDataF E, MyDataF P) {
    MyDataF EDivP;
    EDivP = fabs(E) / P;
    if (EDivP < 50.0) {
        MyDataF val1 = 4.47e-3 * (EDivP)*(EDivP);
        if (EDivP >= 10.0) {
            return val1;
        } else {
            MyDataF val2 = 4.47 / EDivP;
            return (val1 > val2 ? val1 : val2);
        }
    } else {
        MyDataF temp = sqrt(EDivP);
        return (EDivP <= 90 ? 1.58 * temp : 142 / temp);
    }

}
// Calculate We by Nikonov formula

MyDataF We_Nikonov(MyDataF E, MyDataF P) {
    return -0.0382 * E - 2.9e5 * E / P;
}
// Niu_a

MyDataF Niu_a_Nikonov(MyDataF E, MyDataF P) {
    return Eta_Nikonov(E, P) * fabs(We_Nikonov(E, P));
}
// Niu_i

MyDataF Niu_i_Nikonov(MyDataF E, MyDataF P) {
    return Alpha_Nikonov(E, P) * fabs(We_Nikonov(E, P));
}

//////////////////////////////////////////////////////////////////////////
//Morrow and Lowke formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Morrow and Lowke formula

MyDataF Alpha_MorrowAndLowke(MyDataF E, MyDataF N) {
    MyDataF edn = fabs(E) / N / 1e-15;

    if (edn > 1.5) {
        return 2e-16 * N * exp(-7.248 / edn);
    } else {
        return 6.619e-17 * N * exp(-5.593 / edn);
    }
}
// Calculate Eta by Morrow and Lowke formula

MyDataF Eta_MorrowAndLowke(MyDataF E, MyDataF N) {
    MyDataF edn;
    edn = fabs(E) / N / 1e-16;

    if (edn > 0.12) {
        if (edn < 10.50)
            return N * ((6.089e-20 * edn - 2.893e-19) + N * 4.47778e-59 * pow(edn * 1e-16, -1.2749));
        else
            return N * ((8.889e-21 * edn + 2.567e-19) + N * 4.47778e-59 * pow(edn * 1e-16, -1.2749));
    } else {
        if (edn < 0)
            return 0;
        else
            return 106.81;
    }
}
// Calculate We by Morrow and Lowke formula

MyDataF We_MorrowAndLowke(MyDataF E, MyDataF N) {
    MyDataF edn = fabs(E) / N / 1e-16;

    if (edn > 1.0) {
        if (edn <= 20.00)
            return -sign(E)*(1.03e6 * edn + 1.3e6);
        else
            return -sign(E)*(7.4e5 * edn + 7.1e6);
    } else {
        if (edn <= 0.26)
            return -sign(E)*(6.87e6 * edn + 3.38e4);
        else
            return -sign(E)*(7.2973e5 * edn + 1.63e6);
    }
}
// Niu_a

MyDataF Niu_a_MorrowAndLowke(MyDataF E, MyDataF N) {
    return Eta_MorrowAndLowke(E, N) * fabs(We_MorrowAndLowke(E, N));
}
// Niu_i

MyDataF Niu_i_MorrowAndLowke(MyDataF E, MyDataF N) {
    return Alpha_MorrowAndLowke(E, N) * fabs(We_MorrowAndLowke(E, N));
}

//////////////////////////////////////////////////////////////////////////
//Kang formula
//////////////////////////////////////////////////////////////////////////
// Calculate alpha by Kang formula

MyDataF Alpha_Kang(MyDataF E) {
    return 3.5e3 * exp(-1.65e5 / E);
}
// Calculate Eta by Kang formula

MyDataF Eta_Kang(MyDataF E) {
    return 15.0 * exp(-2.5e4 / E);
}
// Calculate We by Kang formula

MyDataF We_Kang(MyDataF E) {
    return -6060.0 * pow(E, 0.75);
}



// Niu_a

MyDataF Niu_a_Kang(MyDataF E) {
    return Eta_Kang(E) * fabs(We_Kang(E));
}
// Niu_i

MyDataF Niu_i_Kang(MyDataF E) {
    return Alpha_Kang(E) * fabs(We_Kang(E));
}

//Calculate Niu_i and Niu_a together

void Niu_Kang(MyDataF *pNiu_i, MyDataF *pNiu_a, MyDataF E) {
    MyDataF We;
    We = fabs(We_Kang(E));
    *pNiu_a = Eta_Kang(E) * We;
    *pNiu_i = Alpha_Kang(E) * We;
}

void Niu_Nikonov(MyDataF *pNiu_i, MyDataF *pNiu_a, MyDataF E, MyDataF P) {
    MyDataF we;
    we = fabs(We_Nikonov(E, P));
    *pNiu_a = Eta_Nikonov(E, P) * we;
    *pNiu_i = Alpha_Nikonov(E, P) * we;
}

void Niu_MorrowAndLowke(MyDataF *pNiu_i, MyDataF *pNiu_a, MyDataF E, MyDataF N) {
    MyDataF we;
    we = fabs(We_MorrowAndLowke(E, N));
    *pNiu_a = Eta_MorrowAndLowke(E, N) * we;
    *pNiu_i = Alpha_MorrowAndLowke(E, N) * we;
}


//==================================================================================================

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
//==================================================================================================

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
//==================================================================================================
/*MyDataF ionization(MyDataF Em, MyDataF p)
{
        Em = Em/100.0;
	

        if(Em/p<54) return 0.0;
        else if(Em/p<=120) return p*(5.0+0.19*(Em/p))*1.0e7*exp(-273.8*p/Em);
        else if(Em/p<3000) return p*54.08*1.0e6*sqrt(Em/p)*exp(-359*p/Em);
        else return 0.0;
}*/

//==================================================================================================

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


//==================================================================================================

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

