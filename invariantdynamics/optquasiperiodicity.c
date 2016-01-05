#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <quadmath.h>
#include <string.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"
//2D phase space functions in L1
__float128 f1(__float128 x, __float128 y) { return cosq(x+y)*cosq(x+y);}
__float128 f2(__float128 x, __float128 y) { return sinq(x+y)*sinq(x+y); }
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

void convergence(int rows, int cols, int time, long double aleastx, long double aleasty,
        long double adeltax, long double adeltay, int fnum, char (*m)[cols]) {
    printf("Running Optimized Version of Quasiperiodicity\n");
    unsigned int i, j, t, v;
    
    FILE *f;
    const char name[] = "outputs/text_optquasi_conv_t%u_g%u_xs%.2Lf_ys%.2Lf_xb%.2Lf_yb%.2Lf.txt";
    char fname[100];
    sprintf(fname, name, time,rows, aleastx, aleasty, aleastx+rows*adeltax,aleasty+cols*adeltay);
    f= fopen(fname,"w");
    fprintf(f,"ONLY BOX IS CERTAIN\n\n");
    unsigned int** traj = (unsigned int**) malloc(sizeof(unsigned int*)*time);
    int ntraj=0;
    for(int o=0; o < time; o++) {
        traj[o] = (unsigned int*) malloc(sizeof(unsigned int)*2);
    }
    __float128  wsum=0;
    for(t=0; t<time; t++) {
        wsum += weight((__float128)t/(__float128)time);
    }

    __float128  x,y, xn, yn;
    __float128  first[fnum];
    __float128  second[fnum];
    __float128  diff[fnum];
    __float128  diff_mag;
    char numzeros;
    __float128  wvar;


    __float128 leastx = aleastx;
    __float128 leasty = aleasty;
    __float128 deltax = adeltax;
    __float128 deltay = adeltay;



    for(i=0; i < rows; i++) {
        printf("i is: %u\n", i);
        for(j=0; j < cols; j++) {
            if(m[i][j]==-1) {
                ntraj++;
                x = leastx + j*deltax+0.5*deltax;
                y = leasty + i*deltay+0.5*deltay;
                memset(first, 0, sizeof(__float128)*fnum);
                for( t=0; t<time; t++) {
                    wvar = weight((__float128)t/(__float128)time);

                    for(v=0; v<fnum;v++) {
                        first[v] += (*fvec[v])(x,y)*wvar;
                    }
                    xn = smod(x+y,2.Q*M_PIq);
                    yn = smod(1.4Q*sinq(x+y)+y,2.Q*M_PIq);
                    traj[t][0] = floorq((x-leastx)/deltax);
                    traj[t][1] = floorq((y-leasty)/deltay);
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
                diff_mag=0.Q;
                for(v=0; v<fnum; v++) {
                    diff[v]=(second[v]-first[v])/wsum;
                    diff_mag+=diff[v]*diff[v];
                }
                diff_mag = sqrtq(diff_mag);
                numzeros = (char)(-1.Q*log10q(diff_mag));
                for(t=0;t<time;t++) {
                    m[(int)fmin(fmax(traj[t][1],0),rows-1)][(int)fmin(fmax(traj[t][0],0),cols-1)] = numzeros;
                }
                //m[i][j]=(unsigned char)numzeros; (old way)
            }
        }
    }
    for(int o=0; o< time; o++) {
        free(traj[o]);
    }
    free(traj);
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            fprintf(f, "m[%u][%u] is %d",i,j,m[i][j]);
        }
    }

    printf("%u\n",m[8][8]);
}

int main() {
    int dim = 50;
    char m[dim][dim];
    for(int i=0; i < dim; i++) {
        for(int j=0; j < dim; j++) {
            m[i][j]=-1;
        }
    }
    for(int i=0; i <1000; i++){
        printf("should be 1: %f\n", cos(i/1000.)*cos(i/1000.) + sin(i/1000.)*sin(i/1000.));
    }

    convergence(dim,dim,1000,0.,0.,2*M_PI/(__float128)dim,2*M_PI/(__float128)dim,1,m);
}
