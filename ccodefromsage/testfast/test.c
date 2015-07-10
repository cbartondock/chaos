#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#define max( a, b ) ( ( a > b) ? a : b ) 
#define min( a, b ) ( ( a < b) ? a : b )

double smod(double x, double y) {
    int div = (int)(x/y);
    double result = x - div*y;
    result = (result < 0) ? result + y : result;
    return result;
}

void forward(
        int maxiter,
        int numsamples,
        int cols,
        int rows,
        double leastx,
        double leasty,
        double deltax,
        double deltay,
        unsigned char (*m)[cols])
{


    int i,j;
    unsigned char keep=1;
    register double x, y, xn, yn;
    double t;
    int row, col;
    int s;
    int count;
    int k=0;
#pragma omp parallel
    while(keep && k<maxiter) {
        keep=0;
        k++; 
        count=0;
#pragma omp parallel for private(i,j,x,y,xn,yn,row,col)
        for(j = 0; j< cols; j++) {
            for(i = 0;i < rows; i++) {
                if(! (m[i][j] & 0xf)) {
                    continue;
                }
                count+=1;
                for(int s = 0; s < numsamples; s++) {
                    for(int s2 = 0; s < numsamples; s++) {
                        //x = leastx+j*deltax +deltax* ((double)rand()/(double)RAND_MAX);
                        //y = leasty+i*deltay + deltay * ((double)rand()/(double)RAND_MAX);
                        x = leastx + j*deltax + deltax*((double)s/(double)numsamples);
                        y = leasty + i*deltay +deltay*((double)s2/(double)numsamples);
                        //xn = smod(x+y,6.283185307);
                        //yn = smod(sin(x+y)+y,6.283185307);

                        xn=1.4 - x*x + .3 * y;
                        yn=x;

                        //xn = .87758 * x - .479426 * y;
                        //yn = .479426 * x + .87758 * y;

                        //xn =10*x;
                        //yn =.5*y;
                        //t= .4-6/(1+x*x+y*y);
                        //xn=1+.98*(x*cos(t)-y*sin(t));
                        //yn=.98*(x*sin(t)+y*cos(t));
                        row = (yn-leasty)/deltay;
                        col = (xn-leastx)/deltax;



                        if(row >= 0 && row < rows && col >= 0 && col < cols && (m[row][col]&0xf) == 1) {
                            m[row][col] = (1<<4)|((m[row][col]<<4)>>4);
                        }
                    }
                }
            }
        }
        printf("num hit forward is %u\n",count);
#pragma omp parallel for private(i,j)
        for(j=0; j< cols; j++) {
            for(i=0; i <rows; i++) {
                if( !keep && ((m[i][j]>>4) != (m[i][j]&0xf))) {
                    keep=1;
                }
                m[i][j] = m[i][j]>> 4;
            }
        }
    }
    if(keep==1){
        printf("negative invariance not reached\n");
    } else {
        printf("negative invariance reached after %u iterations\n",k);
    }
    if(max(rows,cols)<=50){
        for(int i=0; i < rows; i++) {
            printf("[");
            for(int j=0; j<cols;j++) {
                printf("%u,",m[i][j]);
            }
            printf("]\n");
        }
    }

}


int main(int argc, char* argv[]) {
    unsigned char m[10][10];
    for(int i=0; i < 10; ++i) {
        for(int j=0; j<10;++j) {
            m[i][j]=1;
        }
    }
    for(int i=0; i < 10; i++) {
        printf("[");
        for(int j=0; j<10;j++) {
            printf("%u,",m[i][j]);
        }
        printf("]\n");
    }
    //forward(1,1.4,.3,1,1000,50,50,0,0,0.12566,0.12566,m);
    forward(1,100000,10,10,-1,-1,.2,.2,m);

}
