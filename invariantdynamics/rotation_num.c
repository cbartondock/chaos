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


void rotnumber(long double xi1, long double yi1, int time) {
    __float128 x = xi1;
    __float128 y = yi1;
    __float128 xn;
    __float128 yn;
    __float128 xav=0.Q;
    __float128 yav=0.Q;
    __float128 cosval=0.Q;
    __float128 dot, mo, mn, theta;
    __float128 rho=0.Q; 
    __float128 wsum=0.Q;
    for(int t=0; t< time; t++) {
        wsum+=weight((__float128)t/(__float128)time);
    }


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

    printf("The rotation number at (xi=%Lf, yi=%Lf) is %.20Lf\n",xi1,yi1,(long double)rho);
}

int main(char* args) {
    rotnumber(M_PIq+0.001,0.,100000);
}
