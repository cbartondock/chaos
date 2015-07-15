#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"

//2D phase space functions in L1
double f1(double x, double y) { return cos(2*M_PI*y); }
double f2(double x, double y) { return cos(2*M_PI*x)*cos(2*M_PI*y); }
double f3(double x, double y) { return sin(4*M_PI*x)*sin(4*M_PI*y); }
double f4(double x, double y) { return sin(6*M_PI*x)*sin(4*M_PI*y); }
double f5(double x, double y) { return sin(4*M_PI*x)*sin(8*M_PI*y); }
double f6(double x, double y) { return sin(8*M_PI*x)*sin(8*M_PI*y); }
double (*fvec[6]) (double x, double y) = {f1,f2,f3,f4,f5,f6};

void partition(int rows, int cols, int time, double leastx, double leasty,
  double deltax, double deltay, int kgrid,  unsigned long (*m)[cols]) {



    int i, j, t, v;


    double x,y, xn, yn;
    double p[2];
    double evecs[rows][cols][6];
    memset(evecs, 0, sizeof(double)*rows*cols*6);
    for(i=0; i < rows; i++) {
        for(j=0; j < cols; j++) {
            if(m[i][j]==1) {
                x = leastx + j*deltax + 0.5*deltax;
                y = leasty + i*deltay + 0.5*deltay;
                for( t=0; t <time; t++) {
                    for( v=0; v<6;v++) {

                        evecs[i][j][v]+= (*fvec[v])(x,y);
                    }
                    xn = smod(x+y,6.283185307);
                    yn = smod(sin(x+y)+y,6.283185307);
                    //p[0]=x;
                    //p[1]=y;
                    //rk4(p,2*M_PI,.2);
                    //xn = smod(p[0],2*M_PI);
                    //yn = p[1];
                    x = xn;
                    y = yn;
                }
                for( v=0; v<6;v++) {
                    evecs[i][j][v]/=(double)time;
                    evecs[i][j][v]+=1;
                    evecs[i][j][v]/=2.;
                }

            }
        }
    }
    double fdelta = 1./(double)kgrid;
    for(i=0; i < rows; i++) {
        for(j=0; j < cols; j++) {
            if(m[i][j]==1) {
                m[i][j]=0;
                for(v=0; v<6; v++) {

                    m[i][j] += floor(evecs[i][j][v]/fdelta) * pow(kgrid,v);
                }
            }
        }
    }

}

int main() {


}
