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

void curve_part(int weighted, double ax, double ay, double bx, double by, int numpoints, int time, double (*points)[2]) {

    int t, v;
    double wsum=0;
    for(t=0; t<time; t++) {
        wsum += weight((double)t/(double)time);
    }

    double x,y, xn, yn;
    double first;
    double wvar;
    int fnum=1;


    for(int p=0; p < numpoints; p++) {
        x = ax + (bx-ax)*((double)p/(double)numpoints);
        y = ay + (by-ay)*((double)p/(double)numpoints);

        first=0;
        for( t=0; t<time; t++) {
            wvar = weight((double)t/(double)time);

            for(v=0; v<fnum;v++) {
                if(weighted) {
                    first += (*fvec[v])(x,y)*wvar;
                } else {
                    first += (*fvec[v])(x,y);
                }   
            }
            xn = smod(x+y,2*M_PI);
            yn = smod(1.4*sin(x+y)+y,2*M_PI);
            x = xn;
            y = yn;
        }
        
        points[p][0] = (double)p/(double)numpoints;
        if(weighted) {
            points[p][1] = first/wsum;
        } else {
            points[p][1] = first/time;
        }
    }

}

int main() {
}
