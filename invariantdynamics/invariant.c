#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "../sparse_matrix_table/smtable.c"
#define max( a, b ) ( ( a > b) ? a : b )
#define min( a, b ) ( ( a < b) ? a : b )

void calc_invariant(
        int maxiterf,
        int maxiterr,
        int numsamples,
        int cols,
        int rows,
        double leastx,
        double leasty,
        double deltax,
        double deltay,
        unsigned char (*m)[cols])
{

//    FILE *f = fopen("test.txt", "w");
//    FILE *f2 = fopen("test2.txt", "w");
    int i,j;
    unsigned char keep=1;
    double x, y, xn, yn;
    double p[2];
    int row, col;
    int k=0;
    int count;
    while(keep && k<maxiterf) {
        keep=0;
        k++;
        count=0;
        for(j = 0; j< cols; j++) {
            for(i = 0;i < rows; i++) {
                if(! (m[i][j] & 0xf)) {
                    continue;
                }
                count++;
                for(int s = 0; s <= numsamples; s++) {
                    for(int s2 = 0; s2 <= numsamples; s2++){
                        x = leastx + j*deltax + deltax * (double)s2/(double)numsamples;
                        y = leasty + i*deltay + deltay * (double)s/(double)numsamples;
                        //xn = 5*x;
                        //yn = .2*y;
                        //xn =1.4-x*x+.3*y;
                        //yn=x;
                        //xn = smod(x+y,6.283185307);
                        //yn = smod(sin(x+y)+y,6.283185307);
                        p[0]=x;
                        p[1]=y;
                        rk4(p,2*M_PI,.2);
                        xn = smod(p[0],2*M_PI);
                        yn = p[1];

                        //xn = .87758 * x - .479426 * y;
                        //yn = .479426 * x + .87758 * y;
                        //t= .4-6/(1+x*x+y*y);
                        //xn=1+.98*(x*cos(t)-y*sin(t));
                        //yn=.98*(x*sin(t)+y*cos(t));
                        //xn= sqrt(x*x*x*x + y*y*y*y + 2*x*x*y*y)*x - 1/2 * x * x - 1/2 * y * y ;
                        //yn= y * sqrt(x*x*x*x + y*y*y*y + 2*x*x*y*y) + x * y;

                        row = (int)floor((yn-leasty)/deltay);
                        col = (int)floor((xn-leastx)/deltax);
                        if(row >= 0 && row < rows && col >= 0 && col < cols && m[row][col]&0xf) {
                            m[row][col] = (1<<4)|m[row][col];
                        }
                    }
                }
            }
        }
        printf("num hit forward: %u\n",count);
        for(j=0; j< cols; j++) {
            for(i=0; i <rows; i++) {
                if( !keep && ((m[i][j]>>4)!=(m[i][j]&0xf))) {
                    keep=1;
                }
                m[i][j] =m[i][j]>> 4;
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
                printf("%u,",m[i][j]);
            }
            printf("]\n");
        }}
    keep=1; k=0;
    while(keep && k<maxiterr) {
        keep=0;
        k++;
        count=0;
        for(j = 0; j< cols; j++) {
            for(i = 0;i < rows; i++) {

                if(!(m[i][j] & 0xf)) {
                    continue;
                }
                count+=1;
                for(int s=0; s <= numsamples; s++) {
                    for(int s2=0; s2<=numsamples;s2++) {
                        x = leastx + j*deltax + deltax * (double)s2/(double)numsamples;
                        y = leasty + i*deltay + deltay * (double)s/(double)numsamples;
                        //xn = 5*x;
                        //yn = .2*y;
                        //xn =1.4-x*x+.3*y;
                        //yn = x;
                        //xn = smod(x+y,6.283185307);
                        //yn = smod(sin(x+y)+y,6.283185307);
                        //xn = .87758 * x - .479426 * y;
                        //yn = .479426 * x + .87758 * y;
                        p[0]=x;
                        p[1]=y;
                        rk4(p,2*M_PI,.2);
                        xn = smod(p[0],2*M_PI);
                        yn = p[1];
                        //t= .4-6/(1+x*x+y*y);
                        //xn=1+.98*(x*cos(t)-y*sin(t));
                        //yn=.98*(x*sin(t)+y*cos(t));
                        //xn= sqrt(x*x*x*x + y*y*y*y + 2*x*x*y*y)*x - 1/2 * x * x - 1/2 * y * y ;
                        //yn= y * sqrt(x*x*x*x + y*y*y*y + 2*x*x*y*y) + x * y;
                        row = (int)floor((yn-leasty)/deltay);
                        col = (int)floor((xn-leastx)/deltax);

                        if(row>=0 && row<rows && col>=0 && col<cols  && m[row][col]&0xf) {
                            m[i][j] = (1<<4)|m[i][j];
                        }
                    }
                }
            }
        }
        printf("num hit reverse: %u\n",count);
        for(j=0; j< cols; j++) {
            for(i=0; i <rows; i++) {
                if( !keep && ((m[i][j]>>4)!=(m[i][j]&0xf))) {
                    keep=1;
                }
                m[i][j] =m[i][j]>> 4;
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
    calc_invariant(1,1,1000,10,10,-1,-1,.2,.2,m);

}
