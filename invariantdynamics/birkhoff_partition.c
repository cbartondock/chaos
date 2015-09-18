#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include "../rk4/rk4.c"
#include "../usefulfunctions/functions.c"

//rapidly spinning harmonic! Should pick up on any and all distortions.
double base(double x, double y) {return 500*(sin(x)*sin(x)+sin(y)*sin(y));} //ellipse parametrization
double h1(double x, double y) {return cos(base(x,y));} //real part
double h2(double x, double y) {return sin(base(x,y));} //imaginary part
double (*hvec[2]) (double x, double y) = {h1, h2};

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
    return exp((1)/(t*(t-1)));
    //return t*(1-t);
}


void partition(int weighted, int rows, int cols, int time, double leastx, double leasty,
        double deltax, double deltay, int kgrid, unsigned long (*m)[cols]) {

    FILE *f;
    
    const char name[] = "outputs/text_birkhoff_fvals_t%u_g%u_xs%.2f_ys%.2f_xb%.2f_yb%.2f_%s_2.txt";
    char fname[100];
    if(weighted) {
        sprintf(fname, name, time,rows, leastx, leasty, leastx+rows*deltax,leasty+cols*deltay,"weighted");
    } else {
        sprintf(fname, name, time,rows, leastx, leasty, leastx+rows*deltax,leasty+cols*deltay,"unweighted");
    }
    f= fopen(fname,"w");
    clock_t begin,end;
    begin = clock();
    int i, j, t, v;
    double x,y, xn, yn;
    double wsum=0;
    double* weights = (double*) malloc(sizeof(double)*time);
    for(t=0; t<time; t++) {
        weights[t]=weight((double)t/(double)time);
        wsum += weights[t];
    }
    double ***evecs = (double***) malloc(sizeof(double**)*rows);
    for(int i=0; i < rows; i++) {
        evecs[i] = (double **) malloc(sizeof(double*)*cols);
        for(int j=0; j < cols; j++) {
            evecs[i][j] = calloc(sizeof(double),2);
        }
    }
    for(i=0; i < rows; i++) {
        for(j=0; j < cols; j++) {
            if(m[i][j]==1) {
                x = leastx + j*deltax + 0.5*deltax;
                y = leasty + i*deltay + 0.5*deltay;
                for( t=0; t <time; t++) {
                    for( v=0; v<2;v++) {
                        if(weighted) {
                            evecs[i][j][v]+= (*hvec[v])(x,y)*weights[t];
                        } else {
                            evecs[i][j][v]+= (*hvec[v])(x,y)*1.;
                        }
                    }
                    xn = smod(x+y,2*M_PI);
                    yn = smod(sin(x+y)+y,2*M_PI);
                    x = xn;
                    y = yn;
                }
                for( v=0; v<2;v++) {
                    if(weighted) {
                        evecs[i][j][v]/=wsum;
                    } else {
                        evecs[i][j][v]/=(double)time;
                    }
                }

            }
        }
    }
    double fdelta = 1./(double)kgrid;
    for(i=0; i < rows; i++) {
        for(j=0; j < cols; j++) {

            fprintf(f,"i: %u, j: %u,x: %.2f, y: %.2f, h1: %.2f, h2: %.2f, h2/h1: %.2f\n", i, j, leastx+j*deltax,leasty+i*deltay,evecs[i][j][0],evecs[i][j][1],evecs[i][j][1]/evecs[i][j][0]);
            if(m[i][j]==1) {
                m[i][j]=0;
                //for(v=0; v<2; v++) {
                printf("evecs[%u][%u][0] = %f\n",i,j,evecs[i][j][0]);
                printf("evecs[%u][%u][1] = %f\n",i,j,evecs[i][j][1]);
                m[i][j] = floor(100000*evecs[i][j][1]/evecs[i][j][0]);//+= floor(evecs[i][j][v]/fdelta) * pow(kgrid,v);

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
