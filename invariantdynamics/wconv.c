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
    long double val;
    long double wvar;
    long double wsum;
    long double* first;
    long double* second;
    long double firstav;
    long double secondav;
    for(int n=0; n < num; n++) {
        x = x0;
        y = y0;
        wsum=0;
        if(weighted) {
            for(int t=0; t<n*time; t++) {
                wsum += weight((long double)t/(long double)(n*time));
            }
        }
        first = (long double *) malloc(sizeof(long double)*n*time);
        for(int t=0; t < n*time; t++) {
            if(weighted) {
                wvar = weight((long double)t/(long double)(n*time));
                val = pow(cos(x+y),1)*wvar;
            } else {val=pow(cos(x+y),1); }

            first[t] =val;
            xn = smod(x + y, 2*M_PI);
            yn = smod(y + 1.4*sin(x+y), 2*M_PI);
            x = xn;
            y = yn;
        }
        firstav = bucketsum(first, n*time);
        if(weighted) { firstav/=wsum; }
        else { firstav/=(n*time); }

        second=(long double *) malloc(sizeof(long double)*n*time);
        for(int t=0; t < n*time; t++) {
            if(weighted) {
                wvar = weight((long double)t/(long double)(n*time));
                val = pow(cos(x+y),1)*wvar;
            } else {val=pow(cos(x+y),1); }

            second[t]=val;
            xn = smod(x + y, 2*M_PI);
            yn = smod(y + 1.4*sin(x+y), 2*M_PI);
            x = xn;
            y = yn;
        }
        secondav = bucketsum(second, n*time);
        if(weighted) { secondav/=wsum; }
        else { secondav/=(n*time); }
        if(n>0) {
            convdata[n-1] = fabs(secondav - firstav);
        }
        free(first);
        free(second);
    }
}

int main() {

}
