#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"
//2D phase space functions in L1
__float128 f1(__float128 x, __float128 y) { return x*x+y*y;}
__float128 f2(__float128 x, __float128 y) { return cosq(x)*cosq(y); }
__float128 f3(__float128 x, __float128 y) { return sinq(4*M_PIq*x)*sinq(4*M_PI*y); }
__float128 f4(__float128 x, __float128 y) { return sinq(6*M_PIq*x)*sinq(4*M_PI*y); }
__float128 f5(__float128 x, __float128 y) { return sinq(4*M_PIq*x)*sinq(8*M_PI*y); }
__float128 f6(__float128 x, __float128 y) { return sinq(8*M_PIq*x)*sinq(8*M_PI*y); }
__float128 (*fvec[6]) (__float128 x, __float128 y) = {f1,f2,f3,f4,f5,f6};

__float128 weight(__float128 t) {
    if(t<=0.Q || t>=1.Q) {
        return 0.Q;
    }
    return expq((1.Q)/(t*(t-1.Q)));
    //return t*(1-t);
}

void drift(__float128 x_i, __float128 y_i, int time, int total) {

    int t, v;
    __float128 wsum=0;
    for(t=0; t<time; t++) {
        wsum += weight((__float128)t/(__float128)time);
    }
    printf("x_i: %.5Lf, y_i: %.5Lf\n",(long double)x_i,(long double)y_i);
    __float128 x = x_i;
    __float128 y = y_i;
    __float128 xn, yn;
    __float128 first;
    printf("Quad Pi is: %.32Lf\n\n\n\n\n", (long double)M_PIq);
    __float128 av, av2;
    __float128 wvar;
    int fnum=1;
    //  for(int k=0; k < total; k++) {
    av=0.Q;
    for(t=0; t< time; t++) {
        printf("x: %.4Lf, y: %.4Lf, av: %.4Lf\n",(long double)x,(long double)y, (long double)av);
        wvar = weight((__float128)t/(__float128)time);
        for(v=0; v < fnum; v++) {
            av+= (*fvec[v])(x,y)*wvar;
        }
        xn = smod(x+y,2.Q*M_PIq);
        yn = smod(1.4Q*sinq(x+y)+y, 2.Q*M_PIq);
        x=xn;
        y=yn;
    }
    av = av/wsum;
    av2=0.Q;
    for(t=0; t< time; t++) {
        wvar = weight((__float128)t/(__float128)time);
        for(v=0; v < fnum; v++) {
            av2+= (*fvec[v])(x,y)*wvar;
        }
        xn = smod(x+y,2.Q*M_PIq);
        yn = smod(1.4Q*sinq(x+y)+y, 2.Q*M_PIq);
        x=xn;
        y=yn;
    }
    av2 = av2/wsum;
    printf("av: %.50Lf, av2: %.50Lf\n",(long double)av,(long double)av2);
    printf("numzeros: %f\n", (double) -1.Q*log10q(sqrtq((av2-av)*(av2-av))));

    //       av=av/wsum;
    //       if(k==0) {
    //           first=av;
    //       }
    //printf("average minus first %.20f\n", av-first);
//}
}

int main() {
    drift(M_PIq, 0.1Q, 20000, 10000);

}
