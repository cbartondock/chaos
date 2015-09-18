#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../usefulfunctions/functions.c"

void exponents(double** points, int numpoints, int time) {
    int i, t;
    long double x,y, xn,yn;
    long double x2,y2, xn2, yn2;
    long double d0 = .00000000001;
    long double d = d0;
    long double l;
    long double lyapunov;
    
    for(i=0; i < numpoints; i++) {
        x = points[i][0];
        y = points[i][1];
        x2 = points[i][0] + d0;
        y2=y;
        lyapunov=0;
        for(t=0; t<time; t++) {
            xn = smod(x + y, 2*M_PI);
            yn = smod(y + 1.4*sin(x+y), 2*M_PI);
            xn2 = smod(x2 +y2, 2*M_PI);
            yn2 = smod(y2+ 1.4*sin(x2+y2), 2*M_PI);
            x = xn;
            y = yn;
            x2 = xn2;
            y2 = yn2;


            d = sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y));
            l = log(d/d0);
            lyapunov+=l;
            
            x2 = x + (d0/d)*(x2-x);
            y2 = y + (d0/d)*(y2-y);
            
        }
        lyapunov/=time;
        points[i][2] = lyapunov;
    }
    for(int i=0; i < numpoints; i++) {
        //printf("points[%u]: (%f, %f, %f, %f)\n",i,points[i][0],points[i][1],points[i][2], points[i][3]);
    }

}

int main() {
}
