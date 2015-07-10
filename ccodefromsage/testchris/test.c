#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#define max( a, b ) ( ( a > b) ? a : b ) 
#define min( a, b ) ( ( a < b) ? a : b )

double smod(double x, double y) {
    int div = (int)(x/y);
    double result = x - div*y;
    result = (result < 0) ? result + y : result;
    return result;
}

void forward(unsigned char mapnum,
        double param1,
        double param2,
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
    double x, y, xn, yn;
    double tx, ty;
    int row, col;

    while(keep && maxiter>0) {
        keep=0;
        maxiter--;
        for(i = 0;i < rows; i++) {
            for(j = 0; j< cols; j++) {
                if(! (m[i][j] & 0xf)) {
                    continue;
                }
                for(int s = 0; s < numsamples; ++s) {
                    x = leastx+j*deltax + deltax * (rand()%500)/500.;
                    y = leasty+i*deltay + deltay * (rand()%500)/500.;
                    xn = /*smod(x+y,6.283185307);*/param1 - x*x + param2 * y;
                    yn = /*smod(sin(x+y)+y,6.283185307);*/x;
                    //double xn =3*x;
                    //double yn =.5*y;
                    row = (yn-leasty+deltay/2.)/deltay;
                    col = (xn-leastx+deltax/2.)/deltax;



                    if(row >= 0 && row < rows && col >= 0 && col < cols) {
                        m[row][col] = (1<<4)|((m[row][col]<<4)>>4);
                    }
                }
            }
        }

        for(int i=0; i <rows; i++) {
            for(int j=0; j< rows; j++) {
                if( !keep && (m[i][j]>>4) ^ ((m[i][j]<<4)>>4)) {
                    keep=1;
                }
                m[i][j] >>= 4;
            }
        }
    }
    if(maxiter==0){
        printf("invariance not reached\n");
    } else {
        printf("invariance reached\n");
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
    forward(1,1.4,.3,1,100000,10,10,-1,-1,.2,.2,m);

}
