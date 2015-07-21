#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"
//2D phase space functions in L1
double f1(double x, double y) { return sin(x)*sin(x); }
double f2(double x, double y) { return sin(x)*sin(x)*sin(x); }
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
    double first[6] = {0,};
    double diff[6] = {0,};
    double diff_mag;
    unsigned char numzeros;
    double evecs[6];
    double wvar;

    double sinval;
    double cosval;
    double originalsin=0;
    double originalcos=0;

    for(i=0; i < rows; i++) {
        printf("i is: %u\n", i);
        for(j=0; j < cols; j++) {
            if(m[i][j]==1) {
                x = leastx + j*deltax + 0.5*deltax;
                y = leasty + i*deltay + 0.5*deltay;
                //printf("x: %f, y: %f\n",x,y);
                //memset(evecs, 0, sizeof(double)*6);
                sinval=0;
                cosval=0;
                for( t=0; t<time; t++) {
                    wvar = weight((double)t/(double)2*time);
                    //printf("First function\n");
                    sinval+= sin(x)*wvar/wsum;
                    cosval+= cos(x)*wvar/wsum;
                    for(v=0; v<fnum;v++) {
                        //evecs[v] += ((*fvec[v])(x,y)*wvar)/wsum;
                        //printf("Next term in sum: %f %f %f\n",evecs[v],wvar,wsum);
                    }
                    //printf("Second function\n");
                    /*for(u=1; u<2; u++) {
                      evecs[u] += ((*fvec[u])(x,y)*wvar)/wsum;
                    //printf("Next term in sum: %f %f %f\n",evecs[u],wvar, wsum);
                    }*/
                    xn = smod(x+y,6.283185307);
                    yn = smod(1.4*sin(x+y)+y,6.283185307);
                    x = xn;
                    y = yn;
                }
                for(v=0; v<fnum; v++) {
                    //first[v]=evecs[v];
                    //printf("First function result: %f\n",first[v]);
                }
                originalsin = sinval;
                originalcos = cosval;
                /*for(u=1; u<2; u++){
                  first2[u]=evecs[u];
                //printf("Second function result: %f\n",first2[u]);
                }*/
                //double sqsum = first[0]+first2[1];
                //printf("sqsum: %f\n",sqsum);
                //return;         
                for( t=0; t <time; t++) {
                    wvar = weight((double)t/(double)time);
                    sinval+=sin(x)*wvar/wsum;
                    cosval+=cos(x)*wvar/wsum;

                    for(v=0; v<fnum;v++) {
                       // evecs[v]+= ((*fvec[v])(x,y)*weight((double)t/(double)time))/wsum;
                    }
                    xn = smod(x+y,6.283185307);
                    yn = smod(1.4*sin(x+y)+y,6.283185307);
                    x = xn;
                    y = yn;
                }
                diff_mag=0;
                for(v=0; v<fnum; v++) {
                    //diff[v]=(evecs[v]-first[v]);
                    //diff_mag+=diff[v]*diff[v];
                }
                diff_mag = sqrt((sinval-originalsin)*(sinval-originalsin));
                diff_mag = sqrt((cosval-originalcos)*(cosval-originalcos));
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

    convergence(dim,dim,1000,0.,0.,2*M_PI/(double)dim,2*M_PI/(double)dim,1,m);
    /*for(int i=0; i < dim; i++) {
      for(int j=0; j < dim; j++) {
      printf("%u ", m[i][j]);
      }
      printf("\n");
      }*/
}
