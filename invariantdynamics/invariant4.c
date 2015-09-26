#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <quadmath.h>
#include "../sparse_matrix_table/smtable.c"
#define max( a, b ) ( ( a > b) ? a : b )
#define min( a, b ) ( ( a < b) ? a : b )



void calc_invariant(
        int maxiterf,
        int maxiterr,
        int numper,
        int cols,
        int rows,
        int chmap,
        int top,
        long double aleastx,
        long double aleasty,
        long double adeltax,
        long double adeltay,
        unsigned char (*m)[cols])
{
    int i,j,c;
    unsigned char keep=1;
    int k=0;
    int count;
    __float128** corners = (__float128**)malloc(sizeof(__float128*)*4); //tl tr br bl
    for(int r=0; r<4; r++) {
        corners[r] = (__float128*)malloc(sizeof(__float128)*2);
    }
    __float128** newcorners = (__float128**)malloc(sizeof(__float128*)*4); //lr oo
    for(int r=0; r<4; r++) {
        newcorners[r] = (__float128*)malloc(sizeof(__float128)*2);
    }

    __float128*** edges = (__float128***) malloc(sizeof(__float128**)*4);
    for(int r=0; r<4; r++) {
        edges[r] = (__float128**)malloc(sizeof(__float128*)*2);
        for(int r2=0; r2<2; r++) {
            edges[r][r2] = (__float128*)malloc(sizeof(__float128)*2);
        }
    }
    
    __float128 leastx = aleastx;
    __float128 leasty = aleasty;
    __float128 deltax = adeltax;
    __float128 deltay = adeltay;

    __float128** edge = (__float128**) malloc(sizeof(__float128*)*2);
    for(int r=0; r<4; r++) {
        edge[r] = (__float128*)malloc(sizeof(__float128)*2);
    }

    __float128* diff = (__float128*) malloc(sizeof(__float128)*2);
    __float128 epsilon=.001Q;
    __float128 tiv, tjv, ceptxraw, ceptyraw;
    __float128 jvminraw,jvmaxraw,ivminraw,ivmaxraw;
    int jvmin,jvmax,ivmin,ivmax,ivmax2,jvmax2;
    int ceptx,cepty,ceptx2,cepty2;
    while(keep && k<maxiterf) {
        keep=0;
        k++;
        count=0;

        for(i = 0;i < rows; i++) {
            for(j = 0; j< cols; j++) {
                if(! (m[i][j] & 0xf)) {
                    continue;
                }
                count++;
                corners[0][0] = leastx + j*deltax; corners[0][1] = leasty + i*deltay; //tl
                corners[1][0] = leastx + (j+1)* deltax; corners[1][1] = leasty + i*deltay; //tr
                corners[2][0] = leastx + (j+1)*deltax; corners[2][1] = leasty + (i+1)*deltay; //br
                corners[3][0] = leastx + (j)*deltax; corners[3][1] = leasty + (i+1)*deltay; //bl
                for(c=0; c<4; c++) {
                    for(int n=0; n < numper; n++) {
                        if(chmap==1){
                            newcorners[c][1] = smod(sinq(corners[c][0]+corners[c][1])+corners[c][1],2.Q*M_PIq);
                            newcorners[c][0] = smod(corners[c][0]+corners[c][1],2.Q*M_PIq);
                        }
                        if(chmap==2){
                            memcpy(&newcorners[c],&corners[c],sizeof(corners[c]));
                            rk4(newcorners[c],2.Q*M_PIq,.001Q);
                        }
                        if(chmap==3){
                            newcorners[c][0] = 1.4Q - corners[c][0]*corners[c][0] + .3*corners[c][1];
                            newcorners[c][1]=corners[c][0];
                        }
                        if(chmap==4){
                            newcorners[c][0] = corners[c][0]/2.;
                            newcorners[c][1] = corners[c][1]*2.;
                        }
                        corners[c][0] = newcorners[c][0];
                        corners[c][1] = newcorners[c][1];
                    }
                }
                memcpy(&edges[0][0], &newcorners[0],sizeof(edges[0][0]));
                memcpy(&edges[0][1], &newcorners[1],sizeof(edges[0][0]));
                memcpy(&edges[1][0], &newcorners[1],sizeof(edges[0][0]));  //   |---->|
                memcpy(&edges[1][1], &newcorners[2],sizeof(edges[0][0]));  //   |     |
                memcpy(&edges[2][0], &newcorners[2],sizeof(edges[0][0]));  //   |<----|
                memcpy(&edges[2][1], &newcorners[3],sizeof(edges[0][0]));  //
                memcpy(&edges[3][0], &newcorners[3],sizeof(edges[0][0]));
                memcpy(&edges[3][1], &newcorners[0],sizeof(edges[0][0]));

                for(int e=0; e<4; e++) {
                    memcpy(&edge,&edges[e],sizeof(edge));
                    diff[0] = edge[1][0]-edge[0][0];
                    diff[1] = edge[1][1]-edge[0][1];
                    ivminraw = (min(edge[0][1],edge[1][1])-leasty)/deltay;
                    ivmaxraw = (max(edge[0][1],edge[1][1])-leasty)/deltay;
                    jvminraw = (min(edge[0][0],edge[1][0])-leastx)/deltax;
                    jvmaxraw = (max(edge[0][0],edge[1][0])-leastx)/deltax;
                    ivmin = ceil(ivminraw-epsilon);
                    ivmax = floor(ivmaxraw+epsilon);
                    jvmin = ceil(jvminraw-epsilon);
                    jvmax = floor(jvmaxraw+epsilon);
                    if(jvmaxraw!=jvminraw) {
                        for(int jv=jvmin; jv<=jvmax; jv++) {
                            tjv = (leastx+jv*deltax-edge[0][0])/(diff[0]);

                            ceptxraw = (edge[0][0] + diff[0]*tjv-leastx)/deltax;
                            ceptyraw = (edge[0][1] + diff[1]*tjv-leasty)/deltay;
                            ceptx = round(ceptxraw);
                            cepty = floor(ceptyraw+epsilon);
                            ceptx2= ceptx-1;
                            cepty2= cepty-1;
                            if(top==2){ //cylindrical topology
                                ceptx = smod(ceptx,cols);
                                ceptx2 = smod(ceptx2,cols);
                            }
                            if(ceptx>=0 &&ceptx<cols && cepty>=0 && cepty<rows && 0xf&m[cepty][ceptx]) {
                                m[cepty][ceptx] = (1<<4)|m[cepty][ceptx];
                            }
                            if(ceptx2>=0 &&ceptx2<cols && cepty>=0 && cepty<rows && 0xf&m[cepty][ceptx2]) {
                                m[cepty][ceptx2] = (1<<4)|m[cepty][ceptx2];
                            }

                            if(smod(ceptyraw,1)<epsilon || smod(ceptyraw,1)>1-epsilon) {
                                if(cepty2>=0 && cepty2<rows && ceptx >=0 && ceptx<cols && 0xf&m[cepty2][ceptx]) {
                                    m[cepty2][ceptx]=(1<<4)|m[cepty2][ceptx];
                                }
                                if(cepty2>=0 && cepty2<rows && ceptx2 >=0 && ceptx2<cols && 0xf&m[cepty2][ceptx2]) {
                                    m[cepty2][ceptx2]=(1<<4)|m[cepty2][ceptx2];
                                }
                            }
                        }
                    }
                    if(ivmaxraw!=ivminraw) {
                        for(int iv=ivmin; iv<=ivmax; iv++) {
                            tiv = (leasty+iv*deltay-edge[0][1])/(diff[1]);

                            ceptxraw = (edge[0][0] + diff[0]*tiv-leastx)/deltax;
                            ceptyraw = (edge[0][1] + diff[1]*tiv-leasty)/deltay;
                            cepty = round(ceptyraw);
                            ceptx = floor(ceptxraw+epsilon);
                            ceptx2= ceptx-1;
                            cepty2= cepty-1;
                            if(top==2){ //cylindrical topology
                                ceptx = smod(ceptx,cols);
                                ceptx2 = smod(ceptx2,cols);
                            }
                            if(ceptx>=0 &&ceptx<cols && cepty >=0 && cepty<rows && 0xf&m[cepty][ceptx]) {
                                m[cepty][ceptx] = (1<<4)|m[cepty][ceptx];
                            }
                            if(ceptx>=0 &&ceptx<cols && cepty2 >=0 && cepty2<rows && 0xf&m[cepty2][ceptx]) {
                                m[cepty2][ceptx] = (1<<4)|m[cepty2][ceptx];
                            }
                            if(smod(ceptxraw,1)<epsilon || smod(ceptxraw,1)>1-epsilon) {
                                if(cepty>=0 && cepty<rows && ceptx2>=0 && ceptx2<cols && 0xf&m[cepty][ceptx2]) {
                                    m[cepty][ceptx2] = (1<<4)|m[cepty][ceptx2];
                                }
                                if(cepty2>=0 && cepty2<rows && ceptx2>=0 && ceptx2<cols && 0xf&m[cepty2][ceptx2])          {
                                    m[cepty2][ceptx2] = (1<<4)|m[cepty2][ceptx2];
                                }

                            }
                        }
                    }
                    jvmax2=jvmax-1;
                    ivmax2=ivmax-1;
                    if(top==2) {
                        jvmax = smod(jvmax,cols);
                        jvmax2 = smod(jvmax,cols);
                    }
                    if(jvmax < jvmin && ivmax < ivmin  && ivmax>=0 && ivmax<rows && jvmax>=0 && jvmax <cols && 0xf&m[ivmax][jvmax]) {
                        m[ivmax][jvmax] = (1<<4)|m[ivmax][jvmax];
                    }
                    if(jvmax==jvmin && ivmax < ivmin) {
                        if(ivmax>=0 && ivmax<rows && jvmax>=0 && jvmax < cols && 0xf&m[ivmax][jvmax]) {
                            m[ivmax][jvmax]=(1<<4)|m[ivmax][jvmax];
                        }
                        if(ivmax>=0 && ivmax<rows && jvmax2>=0 && jvmax2 < cols && 0xf&m[ivmax][jvmax2]) {
                            m[ivmax][jvmax2]=(1<<4)|m[ivmax][jvmax2];
                        }

                    }
                    if(ivmax==ivmin && jvmax < jvmin) {
                        if(ivmax>=0 && ivmax<rows && jvmax>=0 && jvmax < cols && 0xf&m[ivmax][jvmax]) {
                            m[ivmax][jvmax]=(1<<4)|m[ivmax][jvmax];
                        }
                        if(ivmax2>=0 && ivmax2<rows && jvmax>=0 && jvmax < cols && 0xf&m[ivmax2][jvmax]) {
                            m[ivmax2][jvmax]=(1<<4)|m[ivmax2][jvmax];
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
                count++;
                corners[0][0] = leastx + j*deltax; corners[0][1] = leasty + i*deltay; //tl
                corners[1][0] = leastx + (j+1)* deltax; corners[1][1] = leasty + i*deltay; //tr
                corners[2][0] = leastx + (j+1)*deltax; corners[2][1] = leasty + (i+1)*deltay; //br
                corners[3][0] = leastx + (j)*deltax; corners[3][1] = leasty + (i+1)*deltay; //bl
                for(c=0; c<4; c++) {
                    for(int n=0; n < numper; n++) {
                        if(chmap==1){
                            newcorners[c][1] = smod(sinq(corners[c][0]+corners[c][1])+corners[c][1],2.Q*M_PIq);
                            newcorners[c][0] = smod(corners[c][0]+corners[c][1],2.Q*M_PIq);
                        }
                        if(chmap==2){
                            memcpy(&newcorners[c],&corners[c],sizeof(corners[c]));
                            rk4(newcorners[c],2.Q*M_PIq,.001Q);
                        }
                        if(chmap==3){
                            newcorners[c][0] = 1.4Q - corners[c][0]*corners[c][0] + .3*corners[c][1];
                            newcorners[c][1]=corners[c][0];
                        }
                        if(chmap==4){
                            newcorners[c][0] = corners[c][0]/2.;
                            newcorners[c][1] = corners[c][1]*2.;
                        }
                        corners[c][0] = newcorners[c][0];
                        corners[c][1] = newcorners[c][1];
                    }
                }

                memcpy(&edges[0][0], &newcorners[0],sizeof(edges[0][0]));
                memcpy(&edges[0][1], &newcorners[1],sizeof(edges[0][0]));
                memcpy(&edges[1][0], &newcorners[1],sizeof(edges[0][0]));  //   |---->|
                memcpy(&edges[1][1], &newcorners[2],sizeof(edges[0][0]));  //   |     |
                memcpy(&edges[2][0], &newcorners[2],sizeof(edges[0][0]));  //   |<----|
                memcpy(&edges[2][1], &newcorners[3],sizeof(edges[0][0]));  //
                memcpy(&edges[3][0], &newcorners[3],sizeof(edges[0][0]));
                memcpy(&edges[3][1], &newcorners[0],sizeof(edges[0][0]));

                for(int e=0; e<4; e++) {
                    memcpy(&edge,&edges[e],sizeof(edge));
                    diff[0] = edge[1][0]-edge[0][0];
                    diff[1] = edge[1][1]-edge[0][1];
                    ivminraw = (min(edge[0][1],edge[1][1])-leasty)/deltay;
                    ivmaxraw = (max(edge[0][1],edge[1][1])-leasty)/deltay;
                    jvminraw = (min(edge[0][0],edge[1][0])-leastx)/deltax;
                    jvmaxraw = (max(edge[0][0],edge[1][0])-leastx)/deltax;
                    ivmin = ceil(ivminraw-epsilon);
                    ivmax = floor(ivmaxraw+epsilon);
                    jvmin = ceil(jvminraw-epsilon);
                    jvmax = floor(jvmaxraw+epsilon);
                    if(jvmaxraw!=jvminraw) {
                        for(int jv=jvmin; jv<=jvmax; jv++) {
                            tjv = (leastx+jv*deltax-edge[0][0])/(diff[0]);

                            ceptxraw = (edge[0][0] + diff[0]*tjv-leastx)/deltax;
                            ceptyraw = (edge[0][1] + diff[1]*tjv-leasty)/deltay;
                            ceptx = round(ceptxraw);
                            cepty = floor(ceptyraw+epsilon);
                            ceptx2= ceptx-1;
                            cepty2= cepty-1;
                            if(top==2){ //cylindrical topology
                                ceptx = smod(ceptx,cols);
                                ceptx2 = smod(ceptx2,cols);
                            }
                            if(ceptx>=0 &&ceptx<cols && cepty>=0 && cepty<rows && 0xf&m[cepty][ceptx]) {
                                m[i][j]=(1<<4)|m[i][j];
                            }
                            if(ceptx2>=0 &&ceptx2<cols && cepty>=0 && cepty<rows && 0xf&m[cepty][ceptx2]) {
                                m[i][j]=(1<<4)|m[i][j];

                            }

                            if(smod(ceptyraw,1) <=epsilon || smod(ceptyraw,1)>=1-epsilon) {
                                if(cepty2>=0 && cepty2<rows && ceptx >=0 && ceptx<cols && 0xf&m[cepty2][ceptx]) {
                                    m[i][j]=(1<<4)|m[i][j];

                                }
                                if(cepty2>=0 && cepty2<rows && ceptx2 >=0 && ceptx2<cols && 0xf&m[cepty2][ceptx2]) {
                                    m[i][j]=(1<<4)|m[i][j];

                                }
                            }
                        }
                    }
                    if(ivmaxraw!=ivminraw) {
                        for(int iv=ivmin; iv<=ivmax; iv++) {
                            tiv = (leasty+iv*deltay-edge[0][1])/(diff[1]);

                            ceptxraw = (edge[0][0] + diff[0]*tiv-leastx)/deltax;
                            ceptyraw = (edge[0][1] + diff[1]*tiv-leasty)/deltay;
                            cepty = round(ceptyraw);
                            ceptx = floor(ceptxraw+epsilon);
                            ceptx2= ceptx-1;
                            cepty2= cepty-1;
                            if(top==2){ //cylindrical topology
                                ceptx = smod(ceptx,cols);
                                ceptx2 = smod(ceptx2,cols);
                            }
                            if(ceptx>=0 &&ceptx<cols && cepty >=0 && cepty<rows && 0xf&m[cepty][ceptx]) {
                                m[i][j]=(1<<4)|m[i][j];

                            }
                            if(ceptx>=0 &&ceptx<cols && cepty2 >=0 && cepty2<rows && 0xf&m[cepty2][ceptx]) {
                                m[i][j]=(1<<4)|m[i][j];

                            }
                            if(smod(ceptxraw,1)<=epsilon || smod(ceptxraw,1)>=1-epsilon) {
                                if(cepty>=0 && cepty<rows && ceptx2>=0 && ceptx2<cols && 0xf&m[cepty][ceptx2]) {
                                    m[i][j]=(1<<4)|m[i][j];

                                }
                                if(cepty2>=0 && cepty2<rows && ceptx2>=0 && ceptx2<cols && 0xf&m[cepty2][ceptx2])          {
                                    m[i][j]=(1<<4)|m[i][j];

                                }
                            }
                        }
                    }
                    jvmax2=jvmax-1;
                    ivmax2=ivmax-1;
                    if(top==2) {
                        jvmax = smod(jvmax,cols);
                        jvmax2 = smod(jvmax,cols);
                    }
                    if(jvmax < jvmin && ivmax < ivmin  && ivmax>=0 && ivmax<rows && jvmax>=0 && jvmax <cols && 0xf&m[ivmax][jvmax]) {
                        m[i][j] = (1<<4)|m[i][j];
                    }
                    if(jvmax==jvmin && ivmax < ivmin) {
                        if(ivmax>=0 && ivmax<rows && jvmax>=0 && jvmax < cols && 0xf&m[ivmax][jvmax]) {
                            m[i][j]=(1<<4)|m[i][j];
                        }
                        if(ivmax>=0 && ivmax<rows && jvmax2>=0 && jvmax2 < cols && 0xf&m[ivmax][jvmax2]) {
                            m[i][j]=(1<<4)|m[i][j];
                        }

                    }
                    if(ivmax==ivmin && jvmax < jvmin) {
                        if(ivmax>=0 && ivmax<rows && jvmax>=0 && jvmax < cols && 0xf&m[ivmax][jvmax]) {
                            m[i][j]=(1<<4)|m[i][j];
                        }
                        if(ivmax2>=0 && ivmax2<rows && jvmax>=0 && jvmax < cols && 0xf&m[ivmax2][jvmax]) {
                            m[i][j]=(1<<4)|m[i][j];
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
    int dim = 20;
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
    calc_invariant(10,10,1,dim,dim,4,1,-2.5,-2.5,5./dim,5./dim,m);
}
