#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <float.h>
//#include <omp.h>
#include <math.h>
#include </Users/chrisdock/Documents/chaos/ccodefromsage/rk4/rk4.c>
#include </Users/chrisdock/Documents/chaos/ccodefromsage/sparse_matrix_table/smtable.c>
#define max( a, b ) ( ( a > b) ? a : b )
#define min( a, b ) ( ( a < b) ? a : b )

double smod(double x, double y) {
    int div = (int)(x/y);
    double result = x - div*y;
    result = (result < 0) ? result + y : result;
    return result;
}

double fddpx(double x,double y) {
    double p[2];
    p[0]=x; p[1]=y;
    rk4(p,2*M_PI,.2);
    return smod(p[0],2*M_PI);
}
double fddpy(double x,double y) {
    double p[2];
    p[0]=x; p[1]=y;
    rk4(p,2*M_PI,.2);
    return p[1];
}
double fhenonx(double x,double y) {
    return 1.4-x*x+.3*y;
}
double fhenony(double x,double y) {
    return x;
}

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
    double t;
    double p[2];
    int row, col;
    int k=0;
    int count;
    double corners[4][2]; //tl tr bl br
    double newcorners[4][2]; //lr oo
    double tempx;
    double slopes[4]; // lo1 lo2 o1r o2r
    double bslope;
    double tslope;
    int c,c1,c2;
    int g=0;
    int ti,bi;
    //#pragma omp parallel
    printf("intials, deltax: %f, deltay: %f, leastx: %f, leasty: %f\n",deltax,deltay,leastx,leasty);
    while(keep && k<maxiterf) {
        keep=0;
        k++;
        count=0;
        for(j = 0; j< cols; j++) {
            for(i = 0;i < rows; i++) {
                if(! (m[i][j] & 0xf)) {
                    continue;
                }
                //printf("image of (%u,%u) is:\n\n",i,j);
                count++;
                corners[0][0] = leastx + j*deltax; corners[0][1] = leasty + i*deltay; //tl
                corners[1][0] = leastx + (j+1)* deltax; corners[1][1] = leasty + i*deltay; //tr
                corners[2][0] = leastx + j*deltax; corners[2][1] = leasty + (i+1)*deltay; //bl
                corners[3][0] = leastx + (j+1)*deltax; corners[3][1] = leasty + (i+1)*deltay; //br
                /*printf("Corners are:\n");
                printf("topleft (%f,%f)\n",corners[0][0],corners[0][1]);
                printf("topright (%f,%f)\n",corners[1][0],corners[1][1]);
                printf("bottomleft (%f,%f)\n",corners[2][0],corners[2][1]);
                printf("bottomright (%f,%f)\n",corners[3][0],corners[3][1]);*/
                for(c=0; c<4; c++) {
                    tempx = corners[c][0];
                    corners[c][0] = 1.4 - tempx*tempx + .3*corners[c][1];
                    corners[c][1]=tempx;
                    //rk4(corners[c],2*M_PI,.2);
                    //corners[c][0] = smod(corners[c][0],2*M_PI);
                    //corners[c][0] = 2*corners[c][0];
                    //corners[c][1] = .5*corners[c][1];
                }
                /*printf("New Corners are:\n");
                printf("topleft (%f,%f)\n",corners[0][0],corners[0][1]);
                printf("topright (%f,%f)\n",corners[1][0],corners[1][1]);
                printf("bottomleft (%f,%f)\n",corners[2][0],corners[2][1]);
                printf("bottomright (%f,%f)\n",corners[3][0],corners[3][1]);*/
                newcorners[0][0] = DBL_MAX;
                newcorners[1][0] = -DBL_MAX;
                for(c=0; c<4; c++) {
                    if(corners[c][0] < newcorners[0][0]) {
                        c1=c;
                        newcorners[0][0] = corners[c][0];
                        newcorners[0][1] = corners[c][1];
                    }
                    if(corners[c][0] > newcorners[1][0]) {
                        c2=c;
                        newcorners[1][0] = corners[c][0];
                        newcorners[1][1] = corners[c][1];
                    }
                }
                g=0;
                for(c=0; c<4; c++) {
                    if(c!=c1 && c!=c2) {
                        newcorners[2+g][0]=corners[c][0];
                        newcorners[2+g][1]=corners[c][1];
                        g=1;
                    }
                }
                /*printf("Extreme Corners are:\n");
                printf("left (%f,%f)\n",newcorners[0][0],newcorners[0][1]);
                printf("right (%f,%f)\n",newcorners[1][0],newcorners[1][1]);
                printf("other1 (%f,%f)\n",newcorners[2][0],newcorners[2][1]);
                printf("other2 (%f,%f)\n",newcorners[3][0],newcorners[3][1]);*/
                slopes[0] = (newcorners[2][1]-newcorners[0][1])/(newcorners[2][0]-newcorners[0][0]);
                slopes[1] = (newcorners[3][1]-newcorners[0][1])/(newcorners[3][0]-newcorners[0][0]);
                slopes[2] = (newcorners[1][1]-newcorners[2][1])/(newcorners[1][0]-newcorners[2][0]);
                slopes[3] = (newcorners[1][1]-newcorners[3][1])/(newcorners[1][0]-newcorners[3][0]);
                /*printf("Slopes are:\n");
                printf("l->o1: %f\n",slopes[0]);
                printf("l->o2: %f\n",slopes[1]);
                printf("o1->r: %f\n",slopes[2]);
                printf("o2->r: %f\n",slopes[3]);*/
                if(slopes[0]<slopes[1]) {
                    //o1 is above
                    ti=0;
                    bi=1;
                } else {
                    //o2 is above
                    ti=1;
                    bi=0;
                }

                int xvmin = (int)floor((newcorners[0][0]-leastx)/deltax);
                int xvmax = (int)floor((newcorners[1][0]-leastx)/deltax);
                int ttrans = (int)floor((newcorners[ti+2][0]-leastx)/deltax);
                int btrans = (int)floor((newcorners[bi+2][0]-leastx)/deltax);
                //printf("xvmin is %u, xvmax is %u, ttrans is %u, bttrans is %u\n",xvmin,xvmax,ttrans,btrans);
                int xv=xvmin;
                int yvmin = (int)floor((newcorners[0][1]-leasty)/deltay);
                int yvmax = yvmin;
                while(xv<=xvmax) {
                    int yv = min(yvmax,yvmin);
                    //printf("xvmin: %d, xvmax: %d, xv:%d, yvmin: %d, yvmax: %d\n",xvmin,xvmax,xv,yvmin,yvmax);
                    while(yv <= max(yvmin,yvmax)) {
                        if(yv>=0 && yv<rows && xv>=0 && xv<cols  && m[yv][xv]&0xf) {

                            //printf("yv: %u, xv: %u\n",yv,xv);
                            m[yv][xv] = (1<<4)|m[yv][xv];
                        }
                        yv++;
                    }
                    xv++;
                    int h = (xv-xvmin);
                    if(xv<=ttrans) {
                        yvmin = (int) max(floor((newcorners[0][1]+(h-.1)*deltax*slopes[ti] - leasty)/deltay),floor((newcorners[ti+2][1]-leasty)/deltay));
                    } else {
                        yvmin = (int)max(floor((newcorners[ti+2][1]+(h-.1)*deltax*slopes[ti+2] - leasty)/deltay),floor((newcorners[ti+2][1]-leasty)/deltay));
                    }
                    if(xv <= btrans) {
                        yvmax = (int)min(floor((newcorners[0][1]+(h-.1)*deltax*slopes[bi] - leasty)/deltay),floor((newcorners[bi+2][1]-leasty)/deltay));
                    } else {
                        yvmin = (int)min(floor((newcorners[bi+2][1]+(h-.1)*deltax*slopes[bi+2] - leasty)/deltay),floor((newcorners[bi+2][1]-leasty)/deltay));
                    }

                }
            }
        }

        //#pragma omp parallel for private(i,j,x,y,xn,yn,row,col)
        printf("num hit forward: %u\n",count);
        //#pragma omp parallel for private(i,j)
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
        //#pragma omp parallel for private(i,j,x,y,xn,yn,row,col)
        for(j = 0; j< cols; j++) {
            for(i = 0;i < rows; i++) {
                if(! (m[i][j] & 0xf)) {
                    continue;
                }
                //printf("image of (%u,%u) is:\n\n",i,j);
                count++;
                corners[0][0] = leastx + j*deltax; corners[0][1] = leasty + i*deltay; //tl
                corners[1][0] = leastx + (j+1)* deltax; corners[1][1] = leasty + i*deltay; //tr
                corners[2][0] = leastx + j*deltax; corners[2][1] = leasty + (i+1)*deltay; //bl
                corners[3][0] = leastx + (j+1)*deltax; corners[3][1] = leasty + (i+1)*deltay; //br
                /*printf("Corners are:\n");
                printf("topleft (%f,%f)\n",corners[0][0],corners[0][1]);
                printf("topright (%f,%f)\n",corners[1][0],corners[1][1]);
                printf("bottomleft (%f,%f)\n",corners[2][0],corners[2][1]);
                printf("bottomright (%f,%f)\n",corners[3][0],corners[3][1]);*/
                for(c=0; c<4; c++) {
                    tempx = corners[c][0];
                    corners[c][0] = 1.4 - tempx*tempx + .3*corners[c][1];
                    corners[c][1]=tempx;
                    //rk4(corners[c],2*M_PI,.2);
                    //corners[c][0] = smod(corners[c][0],2*M_PI);
                    //corners[c][0] = 2*corners[c][0];
                    //corners[c][1] = .5*corners[c][1];
                }
                /*printf("New Corners are:\n");
                printf("topleft (%f,%f)\n",corners[0][0],corners[0][1]);
                printf("topright (%f,%f)\n",corners[1][0],corners[1][1]);
                printf("bottomleft (%f,%f)\n",corners[2][0],corners[2][1]);
                printf("bottomright (%f,%f)\n",corners[3][0],corners[3][1]);*/
                newcorners[0][0] = DBL_MAX;
                newcorners[1][0] = -DBL_MAX;
                for(c=0; c<4; c++) {
                    if(corners[c][0] < newcorners[0][0]) {
                        c1=c;
                        newcorners[0][0] = corners[c][0];
                        newcorners[0][1] = corners[c][1];
                    }
                    if(corners[c][0] > newcorners[1][0]) {
                        c2=c;
                        newcorners[1][0] = corners[c][0];
                        newcorners[1][1] = corners[c][1];
                    }
                }
                g=0;
                for(c=0; c<4; c++) {
                    if(c!=c1 && c!=c2) {
                        newcorners[2+g][0]=corners[c][0];
                        newcorners[2+g][1]=corners[c][1];
                        g=1;
                    }
                }
                /*printf("Extreme Corners are:\n");
                printf("left (%f,%f)\n",newcorners[0][0],newcorners[0][1]);
                printf("right (%f,%f)\n",newcorners[1][0],newcorners[1][1]);
                printf("other1 (%f,%f)\n",newcorners[2][0],newcorners[2][1]);
                printf("other2 (%f,%f)\n",newcorners[3][0],newcorners[3][1]);*/
                slopes[0] = (newcorners[2][1]-newcorners[0][1])/(newcorners[2][0]-newcorners[0][0]);
                slopes[1] = (newcorners[3][1]-newcorners[0][1])/(newcorners[3][0]-newcorners[0][0]);
                slopes[2] = (newcorners[1][1]-newcorners[2][1])/(newcorners[1][0]-newcorners[2][0]);
                slopes[3] = (newcorners[1][1]-newcorners[3][1])/(newcorners[1][0]-newcorners[3][0]);
                /*printf("Slopes are:\n");
                printf("l->o1: %f\n",slopes[0]);
                printf("l->o2: %f\n",slopes[1]);
                printf("o1->r: %f\n",slopes[2]);
                printf("o2->r: %f\n",slopes[3]);*/
                if(slopes[0]<slopes[1]) {
                    //o1 is above
                    ti=0;
                    bi=1;
                } else {
                    //o2 is above
                    ti=1;
                    bi=0;
                }
                int xvmin = (int)floor((newcorners[0][0]-leastx)/deltax);
                int xvmax = (int)floor((newcorners[1][0]-leastx)/deltax);
                int ttrans = (int)floor((newcorners[ti+2][0]-leastx)/deltax);
                int btrans = (int)floor((newcorners[bi+2][0]-leastx)/deltax);
                //printf("xvmin is %u, xvmax is %u, ttrans is %u, bttrans is %u\n",xvmin,xvmax,ttrans,btrans);
                int xv=xvmin;
                int yvmin = (int)floor((newcorners[0][1]-leasty)/deltay);
                int yvmax = yvmin;
                while(xv<=xvmax) {
                    int yv = min(yvmax,yvmin);
                    while(yv <= max(yvmin,yvmax)) {
                        if(yv>=0 && yv<rows && xv>=0 && xv<cols  && m[yv][xv]&0xf) {

                            //printf("yv: %u, xv: %u\n",yv,xv);
                            m[i][j] = (1<<4)|m[i][j];
                        }
                        yv++;
                    }
                    xv++;
                    int h = (xv-xvmin);
                    if(xv<=ttrans) {
                        yvmin = (int)max(floor((newcorners[0][1]+(h-.1)*deltax*slopes[ti] - leasty)/deltay),floor((newcorners[ti+2][1]-leasty)/deltay));
                    } else {
                        yvmin = (int)max(floor((newcorners[ti+2][1]+(h-.1)*deltax*slopes[ti+2] - leasty)/deltay),floor((newcorners[ti+2][1]-leasty)/deltay));
                    }
                    if(xv <= btrans) {
                        yvmax = (int)min(floor((newcorners[0][1]+(h-.1)*deltax*slopes[bi] - leasty)/deltay),floor((newcorners[bi+2][1]-leasty)/deltay));
                    } else {
                        yvmin = (int)min(floor((newcorners[bi+2][1]+(h-.1)*deltax*slopes[bi+2] - leasty)/deltay),floor((newcorners[bi+2][1]-leasty)/deltay));
                    }
                }
            }
        }
        printf("num hit reverse: %u\n",count);
        //#pragma omp parallel for private(i,j)
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
    int dim = 5;
    unsigned char m[dim][dim];
    for(int i=0; i < dim; ++i) {
        for(int j=0; j<dim;++j) {
            m[i][j]=1;
        }
    }
    for(int i=0; i < dim; i++) {
        printf("[");
        for(int j=0; j<dim;j++) {
            printf("%u,",m[i][j]);
        }
        printf("]\n");
    }
    calc_invariant(30,30,1000,dim,dim,-2.5,-2.5,5./dim,5./dim,m);

    //calc_invariant(30,30,1000,dim,dim,0,-2,6.28/dim,6./dim,m);

}
