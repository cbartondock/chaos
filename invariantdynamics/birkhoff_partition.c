/* Phase space partitioning algorithm using an n dimensional birkhoff average to separate invariant 
 * curves.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"

//rapidly spinning harmonic! Should pick up on any and all distortions.
__float128 base(__float128 x, __float128 y) {return 500.Q*(sinq(x)*sinq(x)+sinq(y)*sinq(y));} //ellipse parametrization
__float128 h1(__float128 x, __float128 y) {return cosq(base(x,y));} //real part
__float128 h2(__float128 x, __float128 y) {return sinq(base(x,y));} //imaginary part
__float128 (*hvec[2]) (__float128 x, __float128 y) = {h1, h2};

//2D phase space functions in L1
double f1(double x, double y) { return cos(2*M_PI*y); }
double f2(double x, double y) { return cos(2*M_PI*x)*cos(2*M_PI*y); }
double f3(double x, double y) { return sin(4*M_PI*x)*sin(4*M_PI*y); }
double f4(double x, double y) { return sin(6*M_PI*x)*sin(4*M_PI*y); }
double f5(double x, double y) { return sin(4*M_PI*x)*sin(8*M_PI*y); }
double f6(double x, double y) { return sin(8*M_PI*x)*sin(8*M_PI*y); }
double (*fvec[6]) (double x, double y) = {f1,f2,f3,f4,f5,f6};


double weight(double t) {
    if(t<=0 || t>=1) { return 0; }
    return expq((1)/(t*(t-1)));
    //return t*(1-t);
}


void partition(int weighted, int rows, int cols, int time, long double aleastx, long double aleasty,
        long double adeltax, long double adeltay, int kgrid, unsigned long (*m)[cols]) {

    FILE *f;
    
    const char name[] = "outputs/text_birkhoff_fvals_t%u_g%u_xs%.2Lf_ys%.2Lf_xb%.2Lf_yb%.2Lf_%s_2.txt";
    char fname[100];
    if(weighted) {
        sprintf(fname, name, time,rows, aleastx, aleasty, aleastx+rows*adeltax,aleasty+cols*adeltay,"weighted");
    } else {
        sprintf(fname, name, time,rows, aleastx, aleasty, aleastx+rows*adeltax,aleasty+cols*adeltay,"unweighted");
    }
    __float128 leastx = aleastx;
    __float128 leasty = aleasty;
    __float128 deltax = adeltax;
    __float128 deltay = adeltay;

    f= fopen(fname,"w");
    clock_t begin,end;
    begin = clock();
    int i, j, t, v;
    __float128 x,y, xn, yn;
    __float128 wsum=0.;
    __float128* weights = (__float128*) malloc(sizeof(__float128)*time);
    for(t=0; t<time; t++) {
        weights[t]=weight((__float128)t/(__float128)time);
        wsum += weights[t];
    }
    __float128 ***evecs = (__float128***) malloc(sizeof(__float128**)*rows);
    for(int i=0; i < rows; i++) {
        evecs[i] = (__float128 **) malloc(sizeof(__float128*)*cols);
        for(int j=0; j < cols; j++) {
            evecs[i][j] = calloc(sizeof(__float128),2);
        }
    }
    for(i=0; i < rows; i++) {
        for(j=0; j < cols; j++) {
            if(m[i][j]==1) {
                x = leastx + j*deltax + 0.5Q*deltax;
                y = leasty + i*deltay + 0.5Q*deltay;
                for( t=0; t <time; t++) {
                    for( v=0; v<2;v++) {
                        if(weighted) {
                            evecs[i][j][v]+= (*hvec[v])(x,y)*weights[t];
                        } else {
                            evecs[i][j][v]+= (*hvec[v])(x,y);
                        }
                    }
                    xn = smod(x+y,2.Q*M_PIq);
                    yn = smod(.5*sin(x+y)+y,2.Q*M_PIq);
                    x = xn;
                    y = yn;
                }
                for( v=0; v<2;v++) {
                    if(weighted) {
                        evecs[i][j][v]/=wsum;
                    } else {
                        evecs[i][j][v]/=(__float128)time;
                    }
                }

            }
        }
    }
    //__float128 fdelta = 1.Q/(__float128)kgrid;
    for(i=0; i < rows; i++) {
        for(j=0; j < cols; j++) {

            fprintf(f,"i: %u, j: %u,x: %.2f, y: %.2f, h1: %.2f, h2: %.2f, h2/h1: %.2f\n", i, j, (double)leastx+j*(double)deltax,(double)leasty+i*(double)deltay,(double)evecs[i][j][0],(double)evecs[i][j][1],(double)evecs[i][j][1]/(double)evecs[i][j][0]);
            if(m[i][j]==1) {
                m[i][j]=0;
                //for(v=0; v<2; v++) {
                printf("evecs[%u][%u][0] = %f\n",i,j,(double)evecs[i][j][0]);
                printf("evecs[%u][%u][1] = %f\n",i,j,(double)evecs[i][j][1]);
                m[i][j] = floor(100000.Q*evecs[i][j][1]/evecs[i][j][0]);//+= floor(evecs[i][j][v]/fdelta) * pow(kgrid,v);

                //}
            }
        }
    }
    fclose(f);
    end=clock();
    double runtime = (double)(end-begin)/(double)CLOCKS_PER_SEC;
    printf("run time: %f\n",runtime);
    /*for(i=0; i <rows; i++) {
        for(j=0; j<cols; j++) {
            printf("%u", m[i][j]);
        }
        printf("\n");
    }*/

}

int main() {

}
