#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"
//2D phase space functions in L1
double f1(double x, double y) { return sin(x);}
double f2(double x, double y) { return cos(x)*cos(y); }
double f3(double x, double y) { return sin(4*M_PI*x)*sin(4*M_PI*y); }
double f4(double x, double y) { return sin(6*M_PI*x)*sin(4*M_PI*y); }
double f5(double x, double y) { return sin(4*M_PI*x)*sin(8*M_PI*y); }
double f6(double x, double y) { return sin(8*M_PI*x)*sin(8*M_PI*y); }
double (*fvec[6]) (double x, double y) = {f1,f2,f3,f4,f5,f6};

double weight(double t) {
    if(t<=0 || t>=1) {
        return 0;
    }
    return exp((1)/(t*(t-1)));
    //return t*(1-t);
}

void convergence(int rows, int cols, int time, double leastx, double leasty,
        double deltax, double deltay, int fnum, unsigned char (*m)[cols]) {

    int i, j, t, v;

    double wsum=0;
    for(t=0; t<time; t++) {
        wsum += weight((double)t/(double)time);
    }

    double x,y, xn, yn;
    double first[fnum];
    double second[fnum];
    double diff[fnum];
    double diff_mag;
    unsigned char numzeros;
    double wvar;
    


    for(i=0; i < rows; i++) {
        printf("i is: %u\n", i);
        for(j=0; j < cols; j++) {
            if(m[i][j]==1) {
                x = leastx + j*deltax + 0.5*deltax;
                y = leasty + i*deltay + 0.5*deltay;

                //printf("x: %f, y: %f\n",x,y);
                memset(first, 0, sizeof(double)*fnum);
                for( t=0; t<time; t++) {
                    wvar = weight((double)t/(double)time);

                    for(v=0; v<fnum;v++) {
                        first[v] += (*fvec[v])(x,y)*wvar;
                    }
                    xn = smod(x+y,6.283185307);
                    yn = smod(1.4*sin(x+y)+y,6.283185307);
                    x = xn;
                    y = yn;
                }
                
                memset(second,0,sizeof(double)*fnum);

                for( t=0; t <time; t++) {
                    wvar = weight((double)t/(double)time);


                    for(v=0; v<fnum;v++) {
                        second[v]+= (*fvec[v])(x,y)*wvar;
                    }
                    xn = smod(x+y,6.283185307);
                    yn = smod(1.4*sin(x+y)+y,6.283185307);
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

    convergence(dim,dim,1000,0.,0.,2*M_PI/(double)dim,2*M_PI/(double)dim,1,m);
    /*for(int i=0; i < dim; i++) {
      for(int j=0; j < dim; j++) {
      printf("%u ", m[i][j]);
      }
      printf("\n");
      }*/
}
