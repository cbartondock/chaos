#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"
//2D phase space functions in L1
double f1(double x, double y) { return sin(x)*sin(x)+sin(y)*sin(y);/*return cos(x+y);*/}
double f2(double x, double y) { return cos(x)*cos(y); }
double f3(double x, double y) { return sin(4*M_PI*x)*sin(4*M_PI*y); }
double f4(double x, double y) { return sin(6*M_PI*x)*sin(4*M_PI*y); }
double f5(double x, double y) { return sin(4*M_PI*x)*sin(8*M_PI*y); }
double f6(double x, double y) { return sin(8*M_PI*x)*sin(8*M_PI*y); }
double (*fvec[6]) (double x, double y) = {f1,f2,f3,f4,f5,f6};

double weight(double t) {
    if(t<=0 || t>=1) {
        return 0;
    }
    return exp((1)/(t*(t-1)));
    //return t*(1-t);
}

void drift(double x_i, double y_i, int time, int total) {

    int t, v;
    double wsum=0;
    for(t=0; t<time; t++) {
        wsum += weight((double)t/(double)time);
    }

    double x = x_i;
    double y = y_i;
    double xn, yn;
    double first;
    double av;
    double wvar;
    int fnum=1;
    for(int k=0; k < total; k++) {
        av=0;
        for(t=0; t< time; t++) {
            wvar = weight((double)t/(double)time);
            for(v=0; v < fnum; v++) {
                av+= (*fvec[v])(x,y)*wvar;
            }
            xn = smod(x+y,2*M_PI);
            yn = smod(1.4*sin(x+y)+y, 2*M_PI);
            x=xn;
            y=yn;
        }
        av=av/wsum;
        if(k==0) {
            first=av;
        }
        printf("average minus first %.20f\n", av-first);
    }
}

int main() {
    drift(M_PI, 0.1, 100000, 10000);
    
}
