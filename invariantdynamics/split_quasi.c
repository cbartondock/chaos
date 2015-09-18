
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

void quasi_search(double** points, int numpoints, int time, int weighted) {
    printf("numpoints: %u\n", numpoints);
    for(int i=0; i < numpoints; i++) {
        //printf("points[%u]: (%f, %f, %f, %f)\n",i,points[i][0],points[i][1],points[i][2], points[i][3]);
    }
    int i, t;
    long double wsum=0;
    if(weighted) {
        for(t=0; t<time; t++) {
            wsum += weight((long double)t/(long double)time);
        }
    }

    long double x,y, xn,yn;
    long double first=0;
    long double wvar;




    for(i=0; i < numpoints; i++) {
        printf("i: %u\n",i);
        x = points[i][0];
        y = points[i][1];

        first=0;
        for(t=0; t<time; t++) {
            if(weighted) {
                wvar = weight((long double)t/(long double)time);
                first += pow(sin(y),4)*wvar;
            } else {first+=pow(sin(y),4); }
            xn = smod(x + y, 2*M_PI);
            yn = smod(y + 1.4*sin(x+y), 2*M_PI);
            x = xn;
            y = yn;

        }
        if(weighted) {
            points[i][2]=first/wsum;
        } else {
            points[i][2]=first/(float)time;
        }
    }
    for(int i=0; i < numpoints; i++) {
        //printf("points[%u]: (%f, %f, %f, %f)\n",i,points[i][0],points[i][1],points[i][2], points[i][3]);
    }

}

int main() {

}
