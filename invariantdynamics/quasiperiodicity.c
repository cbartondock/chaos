#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include <string.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"
//2D phase space functions in L1
long double f1(long double x, long double y) { return cos(x+y);}
long double f2(long double x, long double y) { return cos(x)*cos(y); }
long double f3(long double x, long double y) { return sin(4*M_PI*x)*sin(4*M_PI*y); }
long double f4(long double x, long double y) { return sin(6*M_PI*x)*sin(4*M_PI*y); }
long double f5(long double x, long double y) { return sin(4*M_PI*x)*sin(8*M_PI*y); }
long double f6(long double x, long double y) { return sin(8*M_PI*x)*sin(8*M_PI*y); }
long double (*fvec[6]) (long double x, long double y) = {f1,f2,f3,f4,f5,f6};

long double weight(long double t) {
    if(t<=0 || t>=1) {
        return 0;
    }
    return exp((1)/(t*(t-1)));
    //return t*(1-t);
}

void convergence(int rows, int cols, int time, long double aleastx, long double aleasty,
        long double adeltax, long double adeltay, int fnum, unsigned char (*m)[cols]) {


    FILE *f;

    const char name[] = "outputs/text_quasi_conv_t%u_g%u_xs%.2Lf_ys%.2Lf_xb%.2Lf_yb%.2Lf.txt";
    char fname[100];
    sprintf(fname, name, time,rows, aleastx, aleasty, aleastx+rows*adeltax,aleasty+cols*adeltay);
    f= fopen(fname,"w");
    unsigned int i, j, t, v;

    __float128  wsum=0;
    for(t=0; t<time; t++) {
        wsum += weight((__float128)t/(__float128)time);
    }

    __float128  x,y, xn, yn;
    __float128  first[fnum];
    __float128  second[fnum];
    __float128  diff[fnum];
    __float128  diff_mag;
    unsigned char numzeros;
    __float128  wvar;


    __float128 leastx = aleastx;
    __float128 leasty = aleasty;
    __float128 deltax = adeltax;
    __float128 deltay = adeltay;



    for(i=0; i < rows; i++) {
        printf("i is: %u\n", i);
        for(j=0; j < cols; j++) {
            if(m[i][j]==1) {
                x = leastx + j*deltax + 0.5*deltax;
                y = leasty + i*deltay + 0.5*deltay;
                fprintf(f,"i: %u, j: %u, x: %Lf, y: %Lf, ",i,j,(long double)x,(long double)y);
                memset(first, 0, sizeof(__float128)*fnum);
                for( t=0; t<time; t++) {
                    wvar = weight((__float128)t/(__float128)time);

                    for(v=0; v<fnum;v++) {
                        first[v] += (*fvec[v])(x,y)*wvar;
                    }
                    xn = smod(x+y,2.Q*M_PIq);
                    yn = smod(1.4Q*sinq(x+y)+y,2.Q*M_PIq);
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
                    yn = smod(1.4Q*sinq(x+y)+y,2.Q*M_PIq);
                    x = xn;
                    y = yn;
                }
                diff_mag=0;
                for(v=0; v<fnum; v++) {
                    diff[v]=(second[v]-first[v])/wsum;
                    diff_mag+=diff[v]*diff[v];
                }
                diff_mag = sqrtq(diff_mag);
                numzeros = (int)(-1.Q*log10q(diff_mag));
                fprintf(f, "numzeros: %.2Lf\n", (long double)(-1.Q*log10q(diff_mag)));
                m[i][j]=(unsigned char)numzeros;
            }
        }
    }
    fclose(f);

}

int main() {
    int dim = 50;
    unsigned char m[dim][dim];
    for(int i=0; i < dim; i++) {
        for(int j=0; j < dim; j++) {
            m[i][j]=1;
        }
    }
    for(int i=0; i <1000; i++){
        printf("should be 1: %f\n", cos(i/1000.)*cos(i/1000.) + sin(i/1000.)*sin(i/1000.));
    }

    convergence(dim,dim,1000,0.,0.,2*M_PI/(long double)dim,2*M_PI/(long double)dim,1,m);
}
