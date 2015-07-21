#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void stickiness() {
    double x = .1;
    double y = .1;
    double epsilon = 2*M_PI/20.;
    double xn, yn, tempx, tempy;
    int totaltime=100000;
    int windowtime=5000;
    int windowoverlap=500;
    for(int t=1; t<=totaltime; t++) {
        for(int t2=1; t2<=totaltime; t2++) {
            if(t2%windowtime == windowtime-windowoverlap) {
                tempx=x;
                tempy=y;
            }
            if(t2%windowtime == 0) {

                x=tempx;
                y=tempy;
            }
            xn = smod(x+y,2*M_PI);
            yn = smod(sin(x+y)+y,2*M_PI);
            x=xn;
            y=yn;
        }
        if(t%windowtime == windowtime-windowoverlap) {
            tempx=x;
            tempy=y;
        }
        if(t%windowtime == 0) {

            x=tempx;
            y=tempy;
        }
        xn = smod(x+y,2*M_PI);
        yn = smod(sin(x+y)+y,2*M_PI);
        x=xn;
        y=yn;
    }

}


int main() {

}
