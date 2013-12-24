
#include <stdio.h>
#include <stdlib.h>

typedef double MyDataF;
#include "connectingInterface.h"

#define N 10

double t_0,tauSquare;
int main(int argc,char *argv[]){
   double t,dt,T;
   int i,n,ns;
   if (argc<=1){
      printf("Usage:\n%s n\n",argv[0]);
      return 0;
   }
   T=1.0/110E9;
   n=atoi(argv[1]);
   if(n<10){
      n=100;
   }
   t_0 = 1.2*T;
   tauSquare=T*T/4;
   dt=T/n;
   t=0;
   i=0;
   ns=N*n;
   while(i<ns){
      printf("%5.3e\n",GaussianSource(t));
      t+=dt;
      i++;
   }
   return 0;
}
