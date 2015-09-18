#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

void convergence(int rows, int cols, int time, long double leastx, long double leasty,
        long double deltax, long double deltay, int fnum, unsigned char (*m)[cols]) {

    int i, j, t, v;

    long double wsum=0;
    for(t=0; t<time; t++) {
        wsum += weight((long double)t/(long double)time);
    }

    long double x,y, xn, yn;
    long double first[fnum];
    long double second[fnum];
    long double diff[fnum];
    long double diff_mag;
    unsigned char numzeros;
    long double wvar;
    


    for(i=0; i < rows; i++) {
        printf("i is: %u\n", i);
        for(j=0; j < cols; j++) {
            if(m[i][j]==1) {
                x = leastx + j*deltax + 0.5*deltax;
                y = leasty + i*deltay + 0.5*deltay;

                //printf("x: %f, y: %f\n",x,y);
                memset(first, 0, sizeof(long double)*fnum);
                for( t=0; t<time; t++) {
                    wvar = weight((long double)t/(long double)time);

                    for(v=0; v<fnum;v++) {
                        first[v] += (*fvec[v])(x,y)*wvar;
                    }
                    xn = smod(x+y,2*M_PI);
                    yn = smod(1.4*sin(x+y)+y,2*M_PI);
                    x = xn;
                    y = yn;
                }
                
                memset(second,0,sizeof(long double)*fnum);

                for( t=0; t <time; t++) {
                    wvar = weight((long double)t/(long double)time);


                    for(v=0; v<fnum;v++) {
                        second[v]+= (*fvec[v])(x,y)*wvar;
                    }
                    xn = smod(x+y,2*M_PI);
                    yn = smod(1.4*sin(x+y)+y,2*M_PI);
                    x = xn;
                    y = yn;
                }
                diff_mag=0;
                for(v=0; v<fnum; v++) {
                    diff[v]=(second[v]-first[v])/wsum;
                    diff_mag+=diff[v]*diff[v];
                }
                numzeros = (int)((-1*log(diff_mag))/log(10));
                m[i][j]=(unsigned char)numzeros;
            }
        }
    }

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
    /*for(int i=0; i < dim; i++) {
      for(int j=0; j < dim; j++) {
      printf("%u ", m[i][j]);
      }
      printf("\n");
      }*/
}
