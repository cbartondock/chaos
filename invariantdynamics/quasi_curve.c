/* This program produces the convergence rates of the birkhoff average
 * along a line segment specified by two points. These convergence rates 
 * strongly differentiate the presence of chaos and quasiperiodicity. The 
 * results are processed and plotted by qcrun.py, or stored in files of the
 * form text_quasi_curve_[].txt and plotted by quasi_curve_plotter.py
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include <string.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"
//2D phase space functions in L1
__float128 f1(__float128 x, __float128 y) { return cosq(x+y);}
__float128 f2(__float128 x, __float128 y) { return cosq(x)*cosq(y); }
__float128 f3(__float128 x, __float128 y) { return sinq(4*M_PI*x)*sinq(4*M_PI*y); }
__float128 f4(__float128 x, __float128 y) { return sinq(6*M_PI*x)*sinq(4*M_PI*y); }
__float128 f5(__float128 x, __float128 y) { return sinq(4*M_PI*x)*sinq(8*M_PI*y); }
__float128 f6(__float128 x, __float128 y) { return sinq(8*M_PI*x)*sinq(8*M_PI*y); }
__float128 (*fvec[6]) (__float128 x, __float128 y) = {f1,f2,f3,f4,f5,f6};

__float128 weight(__float128 t) {
    if(t<=0.Q || t>=1.Q) {
        return 0.Q;
    }
    return expq((1.Q)/(t*(t-1.Q)));
}

void curve_convergence(long double ax1, long double ay1, long double bx1, long double by1, long double anonlin, int numpoints, int time, int fnum, int preshift, double (*points)[2]) {

    FILE *f;

    const char name[] = "outputs/text_quasi_curve_t%u_np%u_ax%.2Lf_ay%.2Lf_bx%.2Lf_by%.2Lf.txt";
    char fname[100];
    sprintf(fname, name, time,numpoints,ax1,ay1,bx1,by1);
    f= fopen(fname,"w");

    int t, v;
    __float128 wsum=0;
    for(t=0; t<time; t++) {
        wsum += weight((__float128)t/(__float128)time);
    }

    __float128 ax = ax1;
    __float128 bx = bx1;
    __float128 ay = ay1;
    __float128 by = by1;

    __float128 x,y, xn, yn;
    __float128 first[fnum];
    __float128 second[fnum];
    __float128 diff[fnum];
    __float128 diff_mag;
    __float128 numzeros;
    __float128 wvar;

    __float128 nonlin=anonlin;
    
    for(int p=0; p < numpoints; p++) {
        x = ax + (bx-ax)*((__float128)p/(__float128)numpoints);
        y = ay + (by-ay)*((__float128)p/(__float128)numpoints);
        fprintf(f,"x: %.10Lf, y: %.10Lf, ", (long double)x, (long double)y);
        for(int s=0; s < preshift; s++) {
            for(t=0; t< time; t++) {
                xn = smod(x+y,2.Q*M_PIq);
                yn = smod(nonlin*sinq(x+y)+y,2.Q*M_PIq);
                x = xn;
                y = yn;
            }
        }

        memset(first, 0, sizeof(__float128)*fnum);
        for( t=0; t<time; t++) {
            wvar = weight((__float128)t/(__float128)time);

            for(v=0; v<fnum;v++) {
                first[v] += (*fvec[v])(x,y)*wvar;
            }
            xn = smod(x+y,2.Q*M_PIq);
            yn = smod(nonlin*sinq(x+y)+y,2.Q*M_PIq);
            x = xn;
            y = yn;
        }

        memset(second,0,sizeof(__float128)*fnum);

        for( t=0; t <time; t++) {
            wvar = weight((__float128)t/(__float128)time);


            for(v=0; v<fnum;v++) {
                second[v]+= (*fvec[v])(x,y)*wvar;
            }
            xn = smod(x+y,2.Q*M_PIq);
            yn = smod(nonlin*sinq(x+y)+y,2.Q*M_PIq);
            x = xn;
            y = yn;
        }
        diff_mag=0.Q;
        for(v=0; v<fnum; v++) {
            diff[v]=(second[v]-first[v])/wsum;
            diff_mag+=diff[v]*diff[v];
        }
        diff_mag = sqrtq(diff_mag);
        numzeros = (__float128)(-1.Q*log10q(diff_mag));
        fprintf(f,"numzeros: %.12Lf\n",(long double)numzeros);
        points[p][0] = (double)p/(__float128)numpoints;
        points[p][1] = (double)numzeros;
        printf("numzeros: %.20Lf\n",(long double)numzeros);
    }

}

int main() {
    /*int dim = 50;
      unsigned char m[dim][dim];
      for(int i=0; i < dim; i++) {
      for(int j=0; j < dim; j++) {
      m[i][j]=1;
      }
      }
      for(int i=0; i <1000; i++){
      printf("should be 1: %f\n", cosq(i/1000.)*cosq(i/1000.) + sinq(i/1000.)*sinq(i/1000.));
      }*/

    //convergence(dim,dim,1000,0.,0.,2*M_PI/(__float128)dim,2*M_PI/(__float128)dim,1,m);
    /*for(int i=0; i < dim; i++) {
      for(int j=0; j < dim; j++) {
      printf("%u ", m[i][j]);
      }
      printf("\n");
      }*/
}
