#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../usefulfunctions/functions.c"

long double weight(long double t) {
    if(t<=0 || t>=1) {
        return 0;
    }
    return exp((1)/(t*(t-1)));
    //return t*(1-t);
}


void convergence(double x0, double y0, int time, int num, int weighted, double convdata[time]) {
    long double x,y, xn,yn;
    long double x2,y2, xn2, yn2;
    long double d0 = .00000000001;
    long double d = d0;
    long double val;
    long double wsum;
    long double wvar;
    long double firstav;
    long double secondav;

    for(int n=0; n < num; n++) {
        x = x0;
        y = y0;
        x2 = x0 + d0;
        y2=y;
        wsum=0;
        if(weighted) {
            for(int t=0; t<n*time; t++) {
                wsum += weight((long double)t/(long double)(n*time));
            }
        }
        long double *first = (long double *) malloc(sizeof(long double)*n*time);
        for(int t=0; t<n*time; t++) {
            printf("t is: %u\n",t);

            xn = smod(x + y, 2*M_PI);
            yn = smod(y + 1.4*sin(x+y), 2*M_PI);
            xn2 = smod(x2 +y2, 2*M_PI);
            yn2 = smod(y2+ 1.4*sin(x2+y2), 2*M_PI);
            x = xn;
            y = yn;
            x2 = xn2;
            y2 = yn2;


            d = sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y));
            if(weighted) {
                wvar = weight((long double)t/(long double)(n*time));
                val = log(d/d0)*wvar;
            } else {val = log(d/d0); }
            first[t] = val;
            x2 = x + (d0/d)*(x2-x);
            y2 = y + (d0/d)*(y2-y);
        }
        firstav = bucketsum(first, n*time);
        if(weighted) { firstav/=wsum; }
        else { firstav/=(n*time); }

        long double *second = (long double *) malloc(sizeof(long double)*n*time);
        for(int t=0; t<n*time; t++) {
            printf("t is: %u\n",t);

            xn = smod(x + y, 2*M_PI);
            yn = smod(y + 1.4*sin(x+y), 2*M_PI);
            xn2 = smod(x2 +y2, 2*M_PI);
            yn2 = smod(y2+ 1.4*sin(x2+y2), 2*M_PI);
            x = xn;
            y = yn;
            x2 = xn2;
            y2 = yn2;


            d = sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y));
            if(weighted) {
                wvar = weight((long double)t/(long double)(n*time));
                val = log(d/d0)*wvar;
            } else {val = log(d/d0); }
            second[t] = val;
            x2 = x + (d0/d)*(x2-x);
            y2 = y + (d0/d)*(y2-y);
        }
        secondav = bucketsum(second, n*time);
        if(weighted) { secondav/=wsum; }
        else { secondav/=(n*time); }
        if(n > 1) {
            convdata[n-2] = fabs(secondav-firstav);
        }
    }

}














/*
   void convergence(double x0, double y0, int time, int num, int weighted, double convdata[time]) {
   long double x,y, xn,yn;
   long double x2,y2, xn2, yn2;
   long double d0 = .00000000001;
   long double d = d0;
   long double average=0;
   long double nminusone;
   long double val;
   long double wsum;
   long double wvar;


   for(int n=0; n < num; n++) {
   x = x0;
   y = y0;
   x2 = x0 + d0;
   y2=y;
   wsum=0;
   if(weighted) {
   for(int t=0; t<n*time; t++) {
   wsum += weight((long double)t/(long double)(n*time));
   }
   }

   nminusone=average;
   average=0;
   for(int t=0; t<n*time; t++) {
   printf("t is: %u\n",t);

   xn = smod(x + y, 2*M_PI);
   yn = smod(y + 1.4*sin(x+y), 2*M_PI);
   xn2 = smod(x2 +y2, 2*M_PI);
   yn2 = smod(y2+ 1.4*sin(x2+y2), 2*M_PI);
   x = xn;
   y = yn;
   x2 = xn2;
   y2 = yn2;


   d = sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y));
   if(weighted) {
   wvar = weight((long double)t/(long double)(n*time));
   val = log(d/d0)*wvar;
   } else {val = log(d/d0); }
   average+= val;
   x2 = x + (d0/d)*(x2-x);
   y2 = y + (d0/d)*(y2-y);
   }
   if(weighted) { average/=wsum; }
   else { average/=(n*time); }
   if(n > 1) {
   convdata[n-2] = fabs(average-nminusone);
   }
   }

   }
   */
int main() {
}
