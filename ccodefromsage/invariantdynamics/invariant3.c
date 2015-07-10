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
    printf("initials, deltax: %f, deltay: %f, leastx: %f, leasty: %f\n",deltax,deltay,leastx,leasty);
    while(keep && k<maxiterf) {
        keep=0;
        k++;
        count=0;

        for(i = 0;i < rows; i++) {
            for(j = 0; j< cols; j++) {
                if(! (m[i][j] & 0xf)) {
                    continue;
                }
                printf("image of (%u,%u) is:\n\n",i,j);
                count++;
                corners[0][0] = leastx + j*deltax; corners[0][1] = leasty + i*deltay; //tl
                corners[1][0] = leastx + (j+1)* deltax; corners[1][1] = leasty + i*deltay; //tr
                corners[2][0] = leastx + j*deltax; corners[2][1] = leasty + (i+1)*deltay; //bl
                corners[3][0] = leastx + (j+1)*deltax; corners[3][1] = leasty + (i+1)*deltay; //br
                printf("Corners are:\n");
                printf("topleft (%f,%f)\n",corners[0][0],corners[0][1]);
                printf("topright (%f,%f)\n",corners[1][0],corners[1][1]);
                printf("bottomleft (%f,%f)\n",corners[2][0],corners[2][1]);
                printf("bottomright (%f,%f)\n",corners[3][0],corners[3][1]);
                for(c=0; c<4; c++) {
                    //double x = corners[c][0];
                    //double y = corners[c][1];
                    //corners[c][1] = smod(sin(x+y)+y,2*M_PI);
                    //corners[c][0] = smod(x+y,2*M_PI);
                    tempx = corners[c][0];
                    corners[c][0] = 1.4 - tempx*tempx + .3*corners[c][1];
                    corners[c][1]=tempx;
                    //rk4(corners[c],2*M_PI,.2);
                    //corners[c][0] = 2*corners[c][0];
                    //corners[c][1] = .5*corners[c][1];
                }
                printf("New Corners are:\n");
                printf("topleft (%f,%f)\n",corners[0][0],corners[0][1]);
                printf("topright (%f,%f)\n",corners[1][0],corners[1][1]);
                printf("bottomleft (%f,%f)\n",corners[2][0],corners[2][1]);
                printf("bottomright (%f,%f)\n",corners[3][0],corners[3][1]);
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
                printf("Extreme Corners are:\n");
                printf("left (%f,%f)\n",newcorners[0][0],newcorners[0][1]);
                printf("right (%f,%f)\n",newcorners[1][0],newcorners[1][1]);
                printf("other1 (%f,%f)\n",newcorners[2][0],newcorners[2][1]);
                printf("other2 (%f,%f)\n",newcorners[3][0],newcorners[3][1]);
                slopes[0] = (newcorners[2][1]-newcorners[0][1])/(newcorners[2][0]-newcorners[0][0]);
                slopes[1] = (newcorners[3][1]-newcorners[0][1])/(newcorners[3][0]-newcorners[0][0]);
                slopes[2] = (newcorners[1][1]-newcorners[2][1])/(newcorners[1][0]-newcorners[2][0]);
                slopes[3] = (newcorners[1][1]-newcorners[3][1])/(newcorners[1][0]-newcorners[3][0]);
                printf("Slopes are:\n");
                printf("l->o1: %f\n",slopes[0]);
                printf("l->o2: %f\n",slopes[1]);
                printf("o1->r: %f\n",slopes[2]);
                printf("o2->r: %f\n",slopes[3]);
                if(slopes[0]<slopes[1]) {
                    //o1 is above
                    ti=0;
                    bi=1;
                } else {
                    //o2 is above
                    ti=1;
                    bi=0;
                }
                printf("ti is %u, bi is %u\n",ti,bi);
                int xvmin = (int)floor((newcorners[0][0]-leastx)/deltax);
                int xvmax = (int)floor((newcorners[1][0]-leastx)/deltax);
                int yvmin = (int)floor((newcorners[0][1]-leasty)/deltay);
                int yvmax = yvmin;             
                int ybound1 = (int)floor((newcorners[ti][1]-leasty)/deltay);
                int ybound2 = (int)floor((newcorners[bi][1]-leasty)/deltay);
                int ttrans = (int)floor((newcorners[ti+2][0]-leastx)/deltax);
                int btrans = (int)floor((newcorners[bi+2][0]-leastx)/deltax);
                int xv=xvmin;
                double tslope1 = slopes[ti];
                double tslope2 = slopes[ti+2];
                double bslope1 = slopes[bi];
                double bslope2 = slopes[bi+2];
                double ct = -smod(newcorners[0][1],deltay) + tslope1*smod(newcorners[0][0],deltax); 
                double cb = -smod(newcorners[0][1],deltay) + bslope1*smod(newcorners[0][0],deltax);
                if(ttrans == xvmin) {
                    yvmin = ybound1;
                } else {
                    yvmin=(int)floor((newcorners[0][1]+(deltax-smod(newcorners[0][0],deltax))*tslope1-leasty)/deltay);
                }
                if(btrans == xvmin) {
                    yvmax = ybound2;
                } else {
                    yvmax=(int)floor((newcorners[0][1]+(deltax-smod(newcorners[0][0],deltax))*bslope1-leasty)/deltay);
                }
                printf("xvmin is %d, xvmax is %d, ttrans is %d, btrans is %d\n",xvmin,xvmax,ttrans,btrans);
                printf("tslope1 is %f, tslope2 is %f, bslope1 is %f, bslope2 is %f\n",tslope1,tslope2,bslope1,bslope2);
                while(xv<=xvmax) {
                    int yv = min(yvmin,yvmax);
                    printf("yvmin: %d, yvmax: %d,xv: %d\n",yvmin,yvmax,xv);
                    int xv2 = (int)smod(xv,cols);
                    xv2=xv;
                    while(yv <= max(yvmax,yvmin)) {
                        if(yv>=0 && yv<rows && xv2>=0 && xv2<cols  && m[yv][xv2]&0xf) {
                            printf("yv: %u, xv: %u\n",yv,xv);
                            m[yv][xv2] = (1<<4)|m[yv][xv2];
                        }
                        yv++;
                    }
                    if(xv==ttrans) {
                        ct = -smod(newcorners[ti][1],deltay) + tslope2*smod(newcorners[ti][0],deltax); 

                    }
                    if(xv==btrans) {
                        cb = -smod(newcorners[bi][1],deltay) + bslope2*smod(newcorners[bi][0],deltax);
                    }
                    if(xv>=ttrans) { ct+=tslope2; }
                    else { ct += tslope1; }
                    if(xv>=btrans) { cb+=bslope2; }
                    else { cb+= bslope1; }
                    yvmin += (int)(ct);
                    yvmax += (int)(cb);
                    if(abs(yvmin)>max(abs(ybound1),abs(ybound2))) {
                        yvmin = min(ybound1,ybound2);
                    }
                    if(abs(yvmax)>max(abs(ybound1),abs(ybound2))) {
                        yvmax = max(ybound1,ybound2);
                    }
                    ct = fmod(ct,1);
                    cb = fmod(cb,1);
                    xv++;
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
                    //double x = corners[c][0];
                    //double y = corners[c][1];
                    //corners[c][1] = smod(sin(x+y)+y,2*M_PI);
                    //corners[c][0] = smod(x+y,2*M_PI);
                    tempx = corners[c][0];
                    corners[c][0] = 1.4 - tempx*tempx + .3*corners[c][1];
                    corners[c][1]=tempx;
                    //rk4(corners[c],2*M_PI,.2);
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
                int ybound1 = (int)floor((newcorners[ti][1]-leasty)/deltay);
                int ybound2 = (int)floor((newcorners[bi][1]-leasty)/deltay);
                int yvmin = (int)floor((newcorners[0][1]-leasty)/deltay);
                int yvmax = yvmin;
                int ttrans = (int)floor((newcorners[ti+2][0]-leastx)/deltax);
                int btrans = (int)floor((newcorners[bi+2][0]-leastx)/deltax);
                int xv=xvmin;
                double tslope1 = slopes[ti];
                double tslope2 = slopes[ti+2];
                double bslope1 = slopes[bi];
                double bslope2 = slopes[bi+2];
                double ct = -smod(newcorners[0][1],deltay) + tslope1*smod(newcorners[0][0],deltax); 
                double cb = -smod(newcorners[0][1],deltay) + bslope1*smod(newcorners[0][0],deltax);
                if(ttrans == xvmin) {
                    yvmin = ybound1;
                } else {
                    yvmin=(int)floor((newcorners[0][1]+(deltax-smod(newcorners[0][0],deltax))*tslope1-leasty)/deltay);
                }
                if(btrans == xvmin) {
                    yvmax = ybound2;
                } else {
                    yvmax=(int)floor((newcorners[0][1]+(deltax-smod(newcorners[0][0],deltax))*bslope1-leasty)/deltay);
                }
                while(xv<=xvmax) {
                    int yv = min(yvmin,yvmax);
                    //printf("xvmin: %d, xvmax: %d, xv:%d, yvmin: %d, yvmax: %d\n",xvmin,xvmax,xv,yvmin,yvmax);
                    int xv2 =(int)smod(xv,cols);
                    xv2=xv;
                    while(yv <= max(yvmax,yvmin)) {
                        if(yv>=0 && yv<rows && xv2>=0 && xv2<cols  && m[yv][xv2]&0xf) {
                            //printf("yv: %u, xv: %u\n",yv,xv);
                            m[i][j] = (1<<4)|m[i][j];
                        }
                        yv++;
                    }
                    if(xv==ttrans) {
                        ct = -smod(newcorners[ti][1],deltay) + tslope2*smod(newcorners[ti][0],deltax); 

                    }
                    if(xv==btrans) {
                        cb = -smod(newcorners[bi][1],deltay) + bslope2*smod(newcorners[bi][0],deltax);
                    }
                    if(xv>=ttrans) { ct+=tslope2; }
                    else { ct += tslope1; }
                    if(xv>=btrans) { cb+=bslope2; }
                    else { cb+= bslope1; }
                    yvmin += (int) (ct);
                    yvmax += (int) (cb);
                    //printf("yvmax is: %u, yvmin is: %u, ybound1 is: %u, ybound2 is: %u\n",yvmax,yvmin,ybound1,ybound2);
                    //printf("abs(yvmin) is %u, abs(yvmax) is %u, max(abs) is %u\n", abs(yvmin),abs(yvmax),max(abs(ybound1),abs(ybound2)));
                    if(abs(yvmin)>max(abs(ybound1),abs(ybound2))) {
                        yvmin = min(ybound1,ybound2);
                    }
                    if(abs(yvmax)>max(abs(ybound1),abs(ybound2))) {
                        yvmax = max(ybound1,ybound2);
                    }
                    //printf("yvmax is now: %u\n",yvmax);
                    ct = fmod(ct,1.);
                    cb = fmod(cb,1.);
                    xv++;
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
    calc_invariant(5,5,1000,dim,dim,-2.5,-2.5,5./dim,5./dim,m);

    //calc_invariant(30,30,1000,dim,dim,0,-2,6.28/dim,6./dim,m);

}
