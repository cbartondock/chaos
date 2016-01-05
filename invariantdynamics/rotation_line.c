#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include "../usefulfunctions/functions.c"

__float128 weight(__float128 t) {
    if(t<=0.Q || t>=1.Q) {
        return 0.Q;
    }
    return expq(1.Q/(t*(t-1.Q)));
}


void rotnumber(long double ax1, long double ay1, long double bx1, long double by1, int numpoints, int time) {


    FILE *f;

    const char name[] = "outputs/text_rotation_curve_t%u_np%u_ax%.2Lf_ay%.2Lf_bx%.2Lf_by%.2Lf.txt";
    char fname[100];
    sprintf(fname, name, time,numpoints,ax1,ay1,bx1,by1);
    f= fopen(fname,"w");

    __float128 ax = ax1;
    __float128 ay = ay1;
    __float128 bx = bx1;
    __float128 by = by1; 
    __float128 x = ax;
    __float128 y = ay;
    __float128 xn;
    __float128 yn;

    __float128 wsum=0.Q;
    for(int t=0; t< time; t++) {
        wsum+=weight((__float128)t/(__float128)time);
    }


    for(int p=0; p < numpoints; p++) {
        __float128 xav=0.Q;
        __float128 yav=0.Q;
        __float128 cosval=0.Q;
        __float128 dot, mo, mn, theta;
        __float128 rho=0.Q;
        x = ax + (bx-ax)*((__float128)p/(__float128)numpoints);
        y = ay + (by-ay)*((__float128)p/(__float128)numpoints);

        //printf("The rotation number at (xi=%Lf, yi=%Lf) is ",(long double)x,(long double)y);
        fprintf(f,"x: %.10Lf, y: %.10Lf, ", (long double)x, (long double)y);
        for(int t=0; t < time; t++) {
            xav+=x*weight((__float128)t/(__float128)time);
            yav+=y*weight((__float128)t/(__float128)time);
            xn = smod(x+y,2.Q*M_PIq);
            yn = smod(y+1.4Q*sinq(x+y),2.Q*M_PIq);
            x=xn;
            y=yn;
        }
        xav/=wsum;
        yav/=wsum;

        mn = sqrtq((x-xav)*(x-xav) + (y-yav)*(y-yav));
        for(int t=0; t< time; t++) {
            xn = smod(x+y,2.Q*M_PIq);
            yn = smod(y+1.4Q*sinq(x+y),2.Q*M_PIq);
            dot = (xn-xav)*(x-xav) + (yn-yav)*(y-yav);
            mo = mn;
            mn = sqrtq((xn-xav)*(xn-xav) + (yn-yav)*(yn-yav));
            cosval = dot/(mo*mn);
            theta = acosq(cosval);
            rho += theta*weight((__float128)t/(__float128)time);
            x=xn;
            y=yn;
        }
        rho/=wsum*2.Q*M_PIq;
        //printf("%.20Lf\n", (long double)rho);
        fprintf(f,"rho: %.12Lf\n",(long double)rho);
    }
}

int main(char* args) {
    double p1 = 0;
    double p2 = .75;
    //rotnumber((1-p1)*M_PIq+p1*4.87,0,(1-p2)*M_PIq+p2*4.87,0,1000, 10000);
    rotnumber(3.8 + .03*.28,0,3.8+.03*.32,0,1000,10000);
}
