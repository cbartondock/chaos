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
        int mul,
        int maxiterf,
        int maxiterr,
        int numsamples,
        int cols,
        int rows,
        double leastx,
        double leasty,
        double deltax,
        double deltay,
        uint64_t (*m)[cols],
        uint64_t(*h)[cols*mul])
{


    int i,j;
    unsigned char keep=1;
    register double x, y, xn, yn;
    double t;
    int row, col, hrow, hcol;
    int s, s2,kf;
    int k=0;
    int count;
#pragma omp parallel
    while(keep && k<maxiterf) {
        keep=0;
        k++;
        count=0;
#pragma omp parallel for private(i,j,x,y,xn,yn,row,col)
        for(j = 0; j< cols; j++) {
            for(i = 0;i < rows; i++) {
                if( (((m[i][j]&0x2))>>1) != 1) {
                    continue;
                }
                count++;
                for(s = 0; s <= numsamples; s++) {
                    for(s2 = 0; s2 <=numsamples;s2++) {
                        x = leastx + j*deltax + deltax * ((double)s/(double)(numsamples));
                        y = leasty + i*deltay + deltay * ((double)s2/(double)(numsamples));
                        //xn = 1.9*x;
                        //yn = .5*y;
                        xn =1.4-x*x+.3*y;
                        yn=x;
                        //xn = smod(x+y,6.283185307);
                        //yn = smod(sin(x+y)+y,6.283185307);

                        //xn = cos(1) * x - sin(1) * y;
                        //yn = sin(1) * x + cos(1) * y;
                        //t= .4-6/(1+x*x+y*y);
                        //xn=1+.98*(x*cos(t)-y*sin(t));
                        //yn=.98*(x*sin(t)+y*cos(t));

                        //xn= sqrt(x*x*x*x + y*y*y*y + 2*x*x*y*y)*x - 1/2 * x * x - 1/2 * y * y ;
                        //yn= y * sqrt(x*x*x*x + y*y*y*y + 2*x*x*y*y) + x * y;
                        //xn= y;
                        //yn = -0.2*x+2.75*y-y*y*y;
                        //xn = x*x-y*y+.9*x+.6013*y;
                        //yn = 2*x*y +2*x+.5*y;
                        row = (int)floor((yn-leasty)/deltay);
                        col = (int)floor((xn-leastx)/deltax);
                        hrow = (int)floor((yn - mul*leasty)/deltay);
                        hcol = (int)floor((xn - mul*leastx)/deltax);
                        if(hrow>=0 && hrow < mul*rows && hcol>=0 && hcol < mul*cols) {
                            h[hrow][hcol]|=((uint64_t)1<<(64-k));
                        }

                        if(row >= 0 && row < rows && col >= 0 && col < cols) {
                            m[row][col] = 0x1 | m[row][col];
                        }
                    }
                }
            }
        }
        printf("num hit forward: %u\n",count);
#pragma omp parallel for private(i,j)
        for(j=0; j< cols; j++) {
            for(i=0; i <rows; i++) {
                if( !keep && (((m[i][j]&0x2)>>1)!=(m[i][j]&0x1))) {
                    keep=1;
                }
                if((m[i][j]&0x2)>>1 ==1) {
                    m[i][j]|=((uint64_t)1<<(64-k));
                }
                m[i][j]=((m[i][j]>>2)<<2)|((m[i][j]<<1)&0x3);

            }
        }

    }
    if(keep==1){
        printf("negative invariance not reached\n");
    } else {
        printf("negative invariance reached after %u iterations\n",k);
    }
    if(max(rows,cols)<=50){
        for(int i=0; i < cols; i++) {
            printf("[");
            for(int j=0; j<rows;j++) {
                printf("%u,",(m[i][j]&0x2)>>1);
            }
            printf("]\n");
        }
    }

    kf=k;
    keep=1; k=0;
    while(keep && k<maxiterr) {
        keep=0;
        k++;
        count=0;
#pragma omp parallel for private(i,j,x,y,xn,yn,row,col)
        for(j = 0; j< cols; j++) {
            for(i = 0;i < rows; i++) {

                if( (((m[i][j]&0x2))>>1) != 1) {
                    continue;
                }
                count++;
                for(s=0; s <numsamples; s++) {
                    for(s2=0; s2 <=numsamples; s2++) {
                        x = leastx + j*deltax + deltax * ((double)rand()/(double)RAND_MAX);
                        y = leasty + i*deltay + deltay * ((double)rand()/(double)RAND_MAX);
                        //xn = 1.9*x;
                        //yn = .5*y;
                        xn =1.4-x*x+.3*y;
                        yn=x;
                        //xn = smod(x+y,6.283185307);
                        //yn = smod(sin(x+y)+y,6.283185307);
                        //xn = cos(1) * x - sin(1) * y;
                        //yn = sin(1) * x + cos(1) * y;
                        //t= .4-6/(1+x*x+y*y);
                        //xn=1+.98*(x*cos(t)-y*sin(t));
                        //yn=.98*(x*sin(t)+y*cos(t));
                        //xn= sqrt(x*x*x*x + y*y*y*y + 2*x*x*y*y)*x - 1/2 * x * x - 1/2 * y * y ;
                        //yn= y * sqrt(x*x*x*x + y*y*y*y + 2*x*x*y*y) + x * y;
                        //xn= y;
                        //yn = -0.2*x+2.75*y-y*y*y;
                        //xn = x*x-y*y+.9*x+.6013*y;
                        //yn = 2*x*y +2*x+.5*y;
                        row = (int)floor((yn-leasty)/deltay);
                        col = (int)floor((xn-leastx)/deltay);
                        if(row>=0 && row<rows && col>=0 && col<cols && (((m[row][col]&0x2))>>1) == 1 ) {
                            m[i][j] = 0x1 | m[i][j];
                            continue;
                        }
                    }
                }
            }
        }
        printf("num hit reverse: %u\n",count);


#pragma omp parallel for private(i,j)
        for(j=0; j< cols; j++) {
            for(i=0; i <rows; i++) {
                if( !keep && (((m[i][j]&0x2)>>1)!=(m[i][j]&0x1)) ) {
                    keep=1;
                }
                if((m[i][j]&0x2)>>1 ==1) {
                    m[i][j]|=((uint64_t)1<<(63-kf-k));
                }
                m[i][j] = ((m[i][j]>>2)<<2)|((m[i][j]<<1)&0x3);

            }
        }
    }

    if(keep==1){
        printf("positive invariance not reached\n");
    } else {
        printf("positive invariance reached after %u iterations\n",k);
    }

    if(max(rows,cols)<=50){
        for(int i=0; i < cols; i++) {
            printf("[");
            for(int j=0; j<rows;j++) {
                printf("%u,",(m[i][j]&0x2)>>1);
            }
            printf("]\n");
        }
    }
    printf("\n\n records \n\n");
    if(max(rows,cols)<=50){
        for(int g=0;g<62;g++) {
            for(int i=0; i < cols; i++) {
                printf("[");
                for(int j=0; j<rows;j++) {
                    printf("%u,",((((m[i][j]>>2)<<2)<<g)>>63));
                }
                printf("]\n");
            }

            printf("\n");
        }
    }

}


int main(int argc, char* argv[]) {
    uint64_t m[10][10];
    uint64_t h[40][40];
    for(int i=0; i < 10; ++i) {
        for(int j=0; j<10;++j) {
            m[i][j]=2;
        }
    }
    for(int i=0; i < 40; ++i) {
        for(int j=0; j<40;++j) {
            h[i][j]=0;
        }
    }


    for(int i=0; i < 10; i++) {
        printf("[");
        for(int j=0; j<10;j++) {
            printf("%u,",(m[i][j]&0x2)>>1);
        }
        printf("]\n");
    }
    forward(1,10,10,1000,10,10,-1,-1,.2,.2,m,h);

}


//xn = smod(x+y,6.283185307);
//yn = smod(sin(x+y)+y,6.283185307);

//xn=param1 - x*x + param2 * y;
//yn=x;

//xn = .87758 * x - .479426 * y;
//yn = .479426 * x + .87758 * y;

//xn =10*x;
//yn =.5*y;
//t= .4-6/(1+x*x+y*y);
//xn=1+.98*(x*cos(t)-y*sin(t));
//yn=.98*(x*sin(t)+y*cos(t));

