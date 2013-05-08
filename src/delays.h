/*
 * Calculate delays of fields for incident
 *
 *
 *
 */

#pragma once

#ifndef __DELAYS_H__
#define __DELAYS_H__

int xs, ys, xe, ye;

MyDataF Source(MyDataF t) {
    //MyDataF t_0;
    //MyDataF tau = 1.0E-11;
    //t_0 =  0.8*tau;//dt*ttstep/10;//||t>1000*dt
    //
    if (t < 0)return 0;
    return cos(omega * t); //-cos(2*M_PI*f*t)*exp(-4*M_PI*pow((t-t_0)/tau,2));//1;//
    //return 0;1.0;
}
//

MyDataF delays(MyDataF px, MyDataF py) {
    return ((px - xs) * cos_phi + (py - ys) * sin_phi);
}

void SetOriginPoint() {
    phi = phi - 2 * M_PI * (floor(phi / (2 * M_PI)));
    if ((phi <= 0.5 * M_PI)) {
        xs = tpis;
        ys = tpjs;
        xe = tpie;
        ye = tpje;
    } else if ((phi <= M_PI)) {
        xs = tpie;
        ys = tpjs;
        xe = tpis;
        ye = tpje;
    } else if ((phi <= 1.5 * M_PI)) {
        xs = tpie;
        ys = tpje;
        xe = tpis;
        ye = tpjs;
    } else {
        xs = tpis;
        ys = tpje;
        xe = tpie;
        ye = tpjs;
    }
    cos_phi = cos(phi);
    sin_phi = sin(phi);
    t_per_cell = dx / c;

    printf("xs = %d\tys = %d\nxe = %d\tye = %d\n", xs, ys, xe, ye);
    printf("sin(phi) = %5.4f\tcos(phi) = %5.4f\n", sin_phi, cos_phi);
    ////system("pause");
}
//Delay arrays
MyDataF **DelayEx, **DelayEy, **DelayEz;
MyDataF **DelayHx, **DelayHy, **DelayHz;
int lenx, leny;
int lenmx, lenmy;

#endif
