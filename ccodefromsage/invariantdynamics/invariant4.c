#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
//#include <omp.h>
#include <math.h>
#include "../sparse_matrix_table/smtable.c"
#define max( a, b ) ( ( a > b) ? a : b )
#define min( a, b ) ( ( a < b) ? a : b )



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
        int cols,
        int rows,
        int chmap,
        int top,
        double leastx,
        double leasty,
        double deltax,
        double deltay,
        unsigned char (*m)[cols])
{

    //    FILE *f = fopen("test.txt", "w");
    //    FILE *f2 = fopen("test2.txt", "w");
    int i,j,c;
    unsigned char keep=1;
    double x, y, xn, yn;
    double t;
    double p[2];
    int row, col;
    int k=0;
    int count;
    double corners[4][2]; //tl tr bl br
    double newcorners[4][2]; //lr oo
    double edges[4][2][2];
    double edge[2][2];
    double diff[2];
    double tiv, tjv, ceptx, cepty;
    int jvmin,jvmax,ivmin,ivmax;
    int iceptx,icepty,iceptx2,icepty2;
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
                printf("\nimage of (%u,%u) is:\n",i,j);
                count++;
                corners[0][0] = leastx + j*deltax; corners[0][1] = leasty + i*deltay; //tl
                corners[1][0] = leastx + (j+1)* deltax; corners[1][1] = leasty + i*deltay; //tr
                corners[2][0] = leastx + j*deltax; corners[2][1] = leasty + (i+1)*deltay; //bl
                corners[3][0] = leastx + (j+1)*deltax; corners[3][1] = leasty + (i+1)*deltay; //br
                printf("corners are: [");
                for(c=0; c<4; c++) {
                    printf("(%f,%f), ", corners[c][1],corners[c][0]);
                }
                printf("]\n");
                for(c=0; c<4; c++) {
                    if(chmap==1){
                        newcorners[c][1] = smod(sin(corners[c][0]+corners[c][1])+corners[c][1],2*M_PI);
                        newcorners[c][0] = smod(corners[c][0]+corners[c][1],2*M_PI);
                    }
                    if(chmap==2){
                        memcpy(&newcorners[c],&corners[c],sizeof(corners[c]));
                        rk4(newcorners[c],2*M_PI,.2);
                    }
                    if(chmap==3){
                        newcorners[c][0] = 1.4 - corners[c][0]*corners[c][0] + .3*corners[c][1];
                        newcorners[c][1]=corners[c][0];
                    }
                    if(chmap==4){
                        newcorners[c][0] = 2*corners[c][0];
                        newcorners[c][1] = .5*corners[c][1];
                    }
                }
                printf("new corners are: [");
                for(c=0; c<4; c++) {
                    printf("(%f,%f), ", newcorners[c][1],newcorners[c][0]);
                }
                printf("]\n");

                memcpy(&edges[0][0], &newcorners[0],sizeof(edges[0][0]));
                memcpy(&edges[0][1], &newcorners[1],sizeof(edges[0][0]));
                memcpy(&edges[1][0], &newcorners[1],sizeof(edges[0][0]));  //   |---->|
                memcpy(&edges[1][1], &newcorners[3],sizeof(edges[0][0]));  //   |     |
                memcpy(&edges[2][0], &newcorners[3],sizeof(edges[0][0]));  //   |<----|
                memcpy(&edges[2][1], &newcorners[2],sizeof(edges[0][0]));  //
                memcpy(&edges[3][0], &newcorners[2],sizeof(edges[0][0]));
                memcpy(&edges[3][1], &newcorners[0],sizeof(edges[0][0]));

                for(int e=0; e<4; e++) {
                    memcpy(&edge,&edges[e],sizeof(edge));
                    diff[0] = edge[1][0]-edge[0][0];
                    diff[1] = edge[1][1]-edge[0][1];
                    printf("edge: (%f,%f)-->(%f,%f), diff: (%f,%f)\n",edge[0][1],edge[0][0],edge[1][1],edge[1][0],diff[1],diff[0]);
                    ivmin = (int)ceil((min(edge[0][1],edge[1][1])-leasty)/deltay);
                    ivmax = (int)floor((max(edge[0][1],edge[1][1])-leasty)/deltay);
                    jvmin = (int)ceil((min(edge[0][0],edge[1][0])-leastx)/deltax);
                    jvmax = (int)floor((max(edge[0][0],edge[1][0])-leastx)/deltax);
                    printf("ivmin: %u, ivmax: %u\n",ivmin,ivmax);
                    for(int iv=ivmin; iv<ivmax; iv++) {
                        tiv = (leasty+iv*deltay-edge[0][1])/(diff[1]);

                        ceptx = edge[0][0] + diff[0]*tiv;
                        cepty = edge[0][1] + diff[1]*tiv;
                        printf("tiv: %f\n",tiv);
                        printf("ceptx: %f, cepty: %f\n",ceptx,cepty);
                        iceptx = (int)floor((ceptx-leastx)/deltax);
                        icepty = (int)round((cepty-leasty)/deltay);
                        iceptx2= iceptx-1;
                        icepty2= icepty-1;
                        if(top==2){
                            //cylinder
                            iceptx = smod(iceptx,cols);
                            iceptx2 = smod(iceptx2,cols);
                        }
                        printf("iv: %d, iceptx: %d, icepty %d\n", iv,iceptx,icepty);
                        if(iceptx>=0 &&iceptx<cols && icepty >=0 && icepty<rows && 0xf&m[icepty][iceptx]) {
                            m[icepty][iceptx] = (1<<4)|m[icepty][iceptx];
                        }
                        if(iceptx>=0 &&iceptx<cols && icepty2 >=0 && icepty2<rows && 0xf&m[icepty2][iceptx]) {
                            m[icepty2][iceptx] = (1<<4)|m[icepty2][iceptx];
                        }
                    }
                    printf("jvmin: %u, jvmax: %u\n",ivmin,ivmax);

                    for(int jv=jvmin; jv<jvmax; jv++) {
                        tjv = (leastx+jv*deltax-edge[0][0])/(diff[0]);

                        ceptx = edge[0][0] + diff[0]*tjv;
                        cepty = edge[0][1] + diff[1]*tjv;

                        iceptx = (int)round((ceptx-leastx)/deltax);
                        icepty = (int)floor((cepty-leasty)/deltay);
                        iceptx2= iceptx-1;
                        icepty2= icepty-1;
                        if(top==2){
                            //cylinder
                            iceptx = smod(iceptx,cols);
                            iceptx2 = smod(iceptx2,cols);
                        }
                        printf("jv: %d, iceptx: %d, icepty %d\n", jv,iceptx,icepty);

                        if(iceptx>=0 &&iceptx<cols && icepty>=0 && icepty<rows && 0xf&m[icepty][iceptx]) {
                            m[icepty][iceptx] = (1<<4)|m[icepty][iceptx];
                        }
                        if(iceptx2>=0 &&iceptx2<cols && icepty>=0 && icepty<rows && 0xf&m[icepty][iceptx2]) {
                            m[icepty][iceptx2] = (1<<4)|m[icepty][iceptx2];
                        }
                    }
                    if(jvmax < jvmin && ivmax < ivmin  && ivmax>=0 && ivmax<rows && jvmax>=0 && jvmax <cols && m[ivmax][jvmax]==1) {
                        printf("crucial case\n");
                        m[ivmax][jvmax] = (1<<4)|m[ivmax][jvmax];
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
        }
    }

    keep=1; k=0;
    while(keep && k<maxiterr) {
        keep=0;
        k++;
        count=0;
        for(i = 0;i < rows; i++) {
            for(j = 0; j< cols; j++) {
                if(! (m[i][j] & 0xf)) {
                    continue;
                }
                printf("\nimage of (%u,%u) is:\n",i,j);
                count++;
                corners[0][0] = leastx + j*deltax; corners[0][1] = leasty + i*deltay; //tl
                corners[1][0] = leastx + (j+1)* deltax; corners[1][1] = leasty + i*deltay; //tr
                corners[2][0] = leastx + j*deltax; corners[2][1] = leasty + (i+1)*deltay; //bl
                corners[3][0] = leastx + (j+1)*deltax; corners[3][1] = leasty + (i+1)*deltay; //br
                printf("corners are: [");
                for(c=0; c<4; c++) {
                    printf("(%f,%f), ", corners[c][1],corners[c][0]);
                }
                printf("]\n");
                for(c=0; c<4; c++) {
                    if(chmap==1){
                        newcorners[c][1] = smod(sin(corners[c][0]+corners[c][1])+corners[c][1],2*M_PI);
                        newcorners[c][0] = smod(corners[c][0]+corners[c][1],2*M_PI);
                    }
                    if(chmap==2){
                        memcpy(&newcorners[c],&corners[c],sizeof(corners[c]));
                        rk4(newcorners[c],2*M_PI,.2);
                    }
                    if(chmap==3){
                        newcorners[c][0] = 1.4 - corners[c][0]*corners[c][0] + .3*corners[c][1];
                        newcorners[c][1]=corners[c][0];
                    }
                    if(chmap==4){
                        newcorners[c][0] = 2*corners[c][0];
                        newcorners[c][1] = .5*corners[c][1];
                    }
                }
                printf("new corners are: [");
                for(c=0; c<4; c++) {
                    printf("(%f,%f), ", newcorners[c][1],newcorners[c][0]);
                }
                printf("]\n");

                memcpy(&edges[0][0], &newcorners[0],sizeof(edges[0][0]));
                memcpy(&edges[0][1], &newcorners[1],sizeof(edges[0][0]));
                memcpy(&edges[1][0], &newcorners[1],sizeof(edges[0][0]));  //   |---->|
                memcpy(&edges[1][1], &newcorners[3],sizeof(edges[0][0]));  //   |     |
                memcpy(&edges[2][0], &newcorners[3],sizeof(edges[0][0]));  //   |<----|
                memcpy(&edges[2][1], &newcorners[2],sizeof(edges[0][0]));  //
                memcpy(&edges[3][0], &newcorners[2],sizeof(edges[0][0]));
                memcpy(&edges[3][1], &newcorners[0],sizeof(edges[0][0]));

                for(int e=0; e<4; e++) {
                    memcpy(&edge,&edges[e],sizeof(edge));
                    diff[0] = edge[1][0]-edge[0][0];
                    diff[1] = edge[1][1]-edge[0][1];
                    printf("edge: (%f,%f)-->(%f,%f), diff: (%f,%f)\n",edge[0][1],edge[0][0],edge[1][1],edge[1][0],diff[1],diff[0]);
                    ivmin = (int)ceil((min(edge[0][1],edge[1][1])-leasty)/deltay);
                    ivmax = (int)floor((max(edge[0][1],edge[1][1])-leasty)/deltay);
                    jvmin = (int)ceil((min(edge[0][0],edge[1][0])-leastx)/deltax);
                    jvmax = (int)floor((max(edge[0][0],edge[1][0])-leastx)/deltax);
                    printf("ivmin: %u, ivmax: %u\n",ivmin,ivmax);
                    for(int iv=ivmin; iv<ivmax; iv++) {
                        tiv = (leasty+iv*deltay-edge[0][1])/(diff[1]);

                        ceptx = edge[0][0] + diff[0]*tiv;
                        cepty = edge[0][1] + diff[1]*tiv;
                        printf("tiv: %f\n",tiv);
                        printf("ceptx: %f, cepty: %f\n",ceptx,cepty);
                        iceptx = (int)floor((ceptx-leastx)/deltax);
                        icepty = (int)round((cepty-leasty)/deltay);
                        iceptx2= iceptx-1;
                        icepty2= icepty-1;
                        if(top==2){
                            //cylinder
                            iceptx = smod(iceptx,cols);
                            iceptx2 = smod(iceptx2,cols);
                        }
                        printf("iv: %d, iceptx: %d, icepty %d\n", iv,iceptx,icepty);
                        if(iceptx>=0 &&iceptx<cols && icepty >=0 && icepty<rows && 0xf&m[icepty][iceptx]) {
                            m[i][j] = (1<<4)|m[i][j];
                        }
                        if(iceptx>=0 &&iceptx<cols && icepty2 >=0 && icepty2<rows && 0xf&m[icepty2][iceptx]) {
                            m[i][j] = (1<<4)|m[i][j];
                        }
                    }
                    printf("jvmin: %u, jvmax: %u\n",ivmin,ivmax);

                    for(int jv=jvmin; jv<jvmax; jv++) {
                        tjv = (leastx+jv*deltax-edge[0][0])/(diff[0]);

                        ceptx = edge[0][0] + diff[0]*tjv;
                        cepty = edge[0][1] + diff[1]*tjv;

                        iceptx = (int)round((ceptx-leastx)/deltax);
                        icepty = (int)floor((cepty-leasty)/deltay);
                        iceptx2= iceptx-1;
                        icepty2= icepty-1;
                        if(top==2){
                            //cylinder
                            iceptx = smod(iceptx,cols);
                            iceptx2 = smod(iceptx2,cols);
                        }
                        printf("jv: %d, iceptx: %d, icepty %d\n", jv,iceptx,icepty);

                        if(iceptx>=0 &&iceptx<cols && icepty>=0 && icepty<rows && 0xf&m[icepty][iceptx]) {
                            m[i][j] = (1<<4)|m[i][j];
                        }
                        if(iceptx2>=0 &&iceptx2<cols && icepty>=0 && icepty<rows && 0xf&m[icepty][iceptx2]) {
                            m[i][j] = (1<<4)|m[i][j];
                        }
                    }
                    if(jvmax <= jvmin && ivmax <= ivmin  && ivmax>=0 && ivmax<rows && jvmax>=0 && jvmax <cols && m[ivmax][jvmax]==1) {
                        printf("crucial case 2\n");
                        m[ivmax][jvmax] = (1<<4)|m[ivmax][jvmax];
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
    calc_invariant(5,5,dim,dim,4,1,-2.5,-2.5,5./dim,5./dim,m);

    //calc_invariant(30,30,dim,dim,0,-2,6.28319/dim,6./dim,m);

}
