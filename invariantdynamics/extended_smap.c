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

void quasi_search(double** points, int numpoints, int time) {
    printf("numpoints: %u\n", numpoints);
    for(int i=0; i < numpoints; i++) {
        //printf("points[%u]: (%f, %f, %f, %f)\n",i,points[i][0],points[i][1],points[i][2], points[i][3]);
    }
    int i, t;
    long double wsum=0;
    for(t=0; t<time; t++) {
        wsum += weight((long double)t/(long double)time);
    }

    long double x,y,z, xn,yn,zn;
    long double first=0;
    long double second=0;
    long double diff;
    long double diff_mag;
    long double numzeros;
    long double wvar;

    long double eps = 30.56;
    long double delta = 1;



    for(i=0; i < numpoints; i++) {
        printf("i: %u\n",i);
        x = points[i][0];
        y = points[i][1];
        z = points[i][2];

        first=0;
        for( t=0; t<time; t++) {
            wvar = weight((long double)t/(long double)time);
            first += cos(x+y)*wvar;
            xn = smod(x + eps*sin(2*M_PI*z) + delta*sin(2*M_PI*y) , 1.);
            yn = smod(y + eps*sin(2*M_PI*z), 1.);
            zn = smod(z + x + eps*sin(2*M_PI*z) + delta*sin(2*M_PI*y), 1.);
            x = xn;
            y = yn;
            z = zn;

        }

        second=0;
        for( t=0; t <time; t++) {
            wvar = weight((long double)t/(long double)time);
            second += cos(x+y)*wvar;
            xn = smod(x + eps*sin(2*M_PI*z) + delta*sin(2*M_PI*y) , 1.);
            yn = smod(y + eps*sin(2*M_PI*z), 1.);
            zn = smod(z + x + eps*sin(2*M_PI*z) + delta*sin(2*M_PI*y), 1.);
            x = xn;
            y = yn;
            z = zn;
        }
        //printf("second: %Le, first: %Le\n",second,first);
        diff=(second-first)/wsum;
        diff_mag=diff*diff;
        //printf("diffmag is: %Le\n",diff_mag);
        numzeros = ((-1*log(diff_mag))/log(10));
        points[i][3]=numzeros;
    }
    for(int i=0; i < numpoints; i++) {
        //printf("points[%u]: (%f, %f, %f, %f)\n",i,points[i][0],points[i][1],points[i][2], points[i][3]);
    }

}

int main() {

}
