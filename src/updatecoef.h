#ifndef _UPDATE_COEF_H_
#define _UPDATE_COEF_H_

void UpdateCoeff() {
    MyDataF gamma, a = 0;
    a = dt * vm / 2;
    gamma = 1 + a;

    InitBeta(&beta, gamma);
    InitCee(beta);
    InitCeh(beta);
    InitCev(beta, alpha);
};

#endif 
