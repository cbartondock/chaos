#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../usefulfunctions/functions.c"
#include "../rk4/rk4.c"
#define max( a, b ) ( ( a > b) ? a : b )
#define min( a, b ) ( ( a < b) ? a : b )

struct adj_element_t {
    int imageindex;
    int domindex;
    int imagenumber;
    struct adj_element_t *next;
    struct adj_element_t *prev;
};
typedef struct adj_element_t adj_element;

typedef struct {
    int domnumber;
    int grid;
    long double leastx;
    long double leasty;
    long double deltax;
    long double deltay;
    adj_element **adjacency_lists;
} sparse_adjacency_matrix;

unsigned char exists_vertex(int index, sparse_adjacency_matrix* matrix) {

    for(int n=0; n < matrix->domnumber; n++) {
        if((matrix->adjacency_lists)[n]->domindex == index) {
            return 1;
        }
    }
    return 0;
}
unsigned char exists_edge(int domindex, int imageindex, sparse_adjacency_matrix* matrix) {
    adj_element* current;
    for(int n=0; n < matrix-> domnumber; n++) {
        if( (matrix->adjacency_lists)[n]->domindex == domindex) {
            current = (matrix->adjacency_lists)[n]->next;
            while(current->imageindex != -1) {
                if(current->imageindex == imageindex) {
                    return 1;
                }
                current = current->next;
            }
            return 0;
        }
    }
    return 0;
}
void add_edge(int domindex, int imageindex, sparse_adjacency_matrix* matrix) {
    adj_element* newelement;

    newelement = (adj_element*) malloc(sizeof(adj_element));
    newelement->domindex = domindex;
    newelement->imageindex = imageindex;
    newelement->imagenumber = -1;
    for(int n=0; n < matrix-> domnumber; n++) {
        if( (matrix->adjacency_lists)[n]->domindex == domindex) {

            newelement->next = (matrix->adjacency_lists)[n];
            newelement->prev = (matrix->adjacency_lists)[n]->prev;
            ((matrix->adjacency_lists)[n]->prev)->next = newelement;
            (matrix->adjacency_lists)[n]->prev = newelement;
            (matrix->adjacency_lists)[n]->imagenumber++;
            break;
        }
    }

}
sparse_adjacency_matrix* initialize_sparse_matrix(int grid,int numper,int chmap,int top, long double aleastx, long double aleasty, long double adeltax, long double adeltay, unsigned char (*m)[grid]) {

    int domnumber=0;
    for(int i=0; i<grid; i++) {
        for(int j=0; j<grid; j++) {
            domnumber+=m[i][j];
        }
    }

    __float128 deltax = adeltax;
    __float128 deltay = adeltay;
    __float128 leastx = aleastx;
    __float128 leasty = aleasty;

    sparse_adjacency_matrix* adjmatrix;
    adjmatrix = (sparse_adjacency_matrix*) malloc(sizeof(sparse_adjacency_matrix));
    adjmatrix->adjacency_lists = (adj_element**) malloc(sizeof(adj_element*)*domnumber);
    adjmatrix->grid = grid;
    adjmatrix->leastx = (long double)leastx;
    adjmatrix->leasty = (long double)leasty;
    adjmatrix->deltax = (long double)deltay;
    adjmatrix->deltay = (long double)deltay;
    adjmatrix->domnumber = domnumber;
    int n = 0;
    adj_element* newelement;
    for(int i=0; i <grid; i++){
        for(int j=0; j <grid; j++){
            if(m[i][j]==1) {
                newelement = (adj_element*) malloc(sizeof(adj_element));
                newelement->domindex = j+i*grid;
                newelement->imageindex = -1;
                newelement->imagenumber = 0;
                newelement->next = newelement;
                newelement->prev = newelement;
                (adjmatrix->adjacency_lists)[n]=newelement;
                n++;
            }
        }
    }
    int i,j,c, domindex, imageindex;
    __float128 corners[4][2];
    __float128 newcorners[4][2]; //tl tr br bl
    __float128 edges[4][2][2];
    __float128 edge[2][2];
    __float128 diff[2];
    __float128 epsilon=.0001;
    __float128 tiv,tjv, ceptxraw, ceptyraw;
    __float128 ivmaxraw,ivminraw,jvmaxraw,jvminraw;
    int ivmax,ivmin,jvmax,jvmin, jvmax2,ivmax2;
    int ceptx, cepty,ceptx2,cepty2;
    for(int b=0; b<domnumber; b++) {
        domindex = (adjmatrix->adjacency_lists)[b]->domindex;
        i=domindex/grid;
        j=domindex%grid;

        corners[0][0] = leastx + j*deltax; corners[0][1] = leasty + i*deltay; //tl
        corners[1][0] = leastx + (j+1)* deltax; corners[1][1] = leasty + i*deltay; //tr
        corners[2][0] = leastx + (j+1)*deltax; corners[2][1] = leasty + (i+1)*deltay; //br
        corners[3][0] = leastx + (j)*deltax; corners[3][1] = leasty + (i+1)*deltay; //bl
        for(int n=0; n<numper; n++) {
            for(c=0; c<4; c++) {
                if(chmap==1){
                    newcorners[c][1] = smod(sinq(corners[c][0]+corners[c][1])+corners[c][1],2*M_PIq);
                    newcorners[c][0] = smod(corners[c][0]+corners[c][1],2*M_PIq);
                }
                if(chmap==2){
                    memcpy(&newcorners[c],&corners[c],sizeof(corners[c]));
                    rk4(newcorners[c],2*M_PIq,.001Q);
                }
                if(chmap==3){
                    newcorners[c][0] = 1.4Q - corners[c][0]*corners[c][0] + .3Q*corners[c][1];
                    newcorners[c][1]=corners[c][0];
                }
                if(chmap==4){
                    newcorners[c][0] = 2*corners[c][0];
                    newcorners[c][1] = .5*corners[c][1];
                }
                corners[c][0]=newcorners[c][0];
                corners[c][1]=newcorners[c][1];
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
            ivmin = ceilq(ivminraw-epsilon);
            ivmax = floorq(ivmaxraw+epsilon);
            jvmin = ceilq(jvminraw-epsilon);
            jvmax = floorq(jvmaxraw+epsilon);
            if(jvmaxraw!=jvminraw) {
                for(int jv=jvmin; jv<=jvmax; jv++) {
                    tjv = (leastx+jv*deltax-edge[0][0])/(diff[0]);

                    ceptxraw = (edge[0][0] + diff[0]*tjv-leastx)/deltax;
                    ceptyraw = (edge[0][1] + diff[1]*tjv-leasty)/deltay;
                    ceptx = roundq(ceptxraw);
                    cepty = floorq(ceptyraw+epsilon);
                    ceptx2= ceptx-1;
                    cepty2= cepty-1;
                    if(top==2){//cylindrical topology
                        ceptx = smod(ceptx,grid);
                        ceptx2 = smod(ceptx2,grid);
                    }
                    imageindex = ceptx + grid*cepty;
                    if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                        add_edge(domindex,imageindex,adjmatrix);
                    }
                    imageindex = ceptx2 + grid*cepty;

                    if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                        add_edge(domindex,imageindex,adjmatrix);
                    }

                    if(smod(ceptyraw,1)<=epsilon || smod(ceptyraw,1)>=1-epsilon) {
                        imageindex = ceptx + grid*cepty2;

                        if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                            add_edge(domindex,imageindex,adjmatrix);
                        }
                        imageindex = ceptx2 + grid*cepty2;

                        if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                            add_edge(domindex,imageindex,adjmatrix);
                        }
                    }
                }
            }
            if(ivmaxraw!=ivminraw) {
                for(int iv=ivmin; iv<=ivmax; iv++) {
                    tiv = (leasty+iv*deltay-edge[0][1])/(diff[1]);

                    ceptxraw = (edge[0][0] + diff[0]*tiv-leastx)/deltax;
                    ceptyraw = (edge[0][1] + diff[1]*tiv-leasty)/deltay;
                    cepty = roundq(ceptyraw);
                    ceptx = floorq(ceptxraw+epsilon);
                    ceptx2= ceptx-1;
                    cepty2= cepty-1;
                    if(top==2){ //cylindrical topology
                        ceptx = smod(ceptx,grid);
                        ceptx2 = smod(ceptx2,grid);
                    }
                    imageindex = ceptx + grid*cepty;
                    if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                        add_edge(domindex,imageindex,adjmatrix);
                    }
                    imageindex = ceptx + grid*cepty2;
                    if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                        add_edge(domindex,imageindex,adjmatrix);
                    }
                    if(smod(ceptxraw,1)<=epsilon || smod(ceptxraw,1)>=1-epsilon) {
                        imageindex = ceptx2 + grid*cepty;
                        if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                            add_edge(domindex,imageindex,adjmatrix);
                        }
                        imageindex = ceptx2 + grid*cepty2;
                        if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                            add_edge(domindex,imageindex,adjmatrix);
                        }

                    }
                }
            }
            jvmax2=jvmax-1;
            ivmax2=ivmax-1;
            if(top==2) {
                jvmax = smod(jvmax,grid);
                jvmax2 = smod(jvmax,grid);
            }
            if(jvmax < jvmin && ivmax < ivmin) {
                imageindex = jvmax + grid*ivmax;
                if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                    add_edge(domindex,imageindex,adjmatrix);
                }
            }
            if(jvmax==jvmin && ivmax < ivmin) {
                imageindex = jvmax+ivmax*grid;
                if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                    add_edge(domindex,imageindex,adjmatrix);
                }
                imageindex=jvmax2+grid*ivmax;
                if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                    add_edge(domindex,imageindex,adjmatrix);
                }

            }
            if(ivmax==ivmin && jvmax < jvmin) {
                imageindex = jvmax + grid*ivmax;
                if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                    add_edge(domindex,imageindex,adjmatrix);
                }
                imageindex = jvmax + grid*ivmax2;
                if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                    add_edge(domindex,imageindex,adjmatrix);
                }
            }
        }
    }
    return adjmatrix;
}

void free_matrix(sparse_adjacency_matrix* matrix) {
    adj_element* current;
    adj_element* next;
    for(int n=0; n<matrix->domnumber; n++) {
        current = (matrix->adjacency_lists)[n]->next;
        while(current->imageindex!=-1) {
            next = current->next;
            free(current);
            current = next;
        }
        free(current);
    }
    free(matrix->adjacency_lists);
    free(matrix);
}


adj_element* get_edges(int domindex, sparse_adjacency_matrix* matrix) {
    for(int n=0; n < matrix-> domnumber; n++) {
        if( (matrix->adjacency_lists)[n]->domindex == domindex) {
            return (matrix->adjacency_lists)[n];
        }
    }
    return 0;
}

int num_edges(int domindex, sparse_adjacency_matrix* matrix) {
    for(int n=0; n < matrix-> domnumber; n++) {
        if( (matrix->adjacency_lists)[n]->domindex == domindex) {
            return (matrix->adjacency_lists)[n]->imagenumber;
        }
    }
    return -1;
}
int get_grid(sparse_adjacency_matrix* matrix) {
    return matrix->grid;
}
int get_domnumber(sparse_adjacency_matrix* matrix) {
    return matrix->domnumber;
}

/*
   int main(char* args) {
   printf("running main\n");
   int grid=2;
   unsigned char  (*m)[grid];
   m[0][0]=1; m[0][1]=1; m[1][0]=1; m[1][1] = 1;
   sparse_adjacency_matrix* sam = initialize_sparse_matrix(grid,2,5,-2,-2,2,2,m);
//sam = quadruple_precision(quadruple_precision(sam));
printf("sam's edges are: [");
adj_element* current;
for(int i=0; i <sam->domnumber; i++) {
current = (sam->adjacency_lists)[i]->next;
while(current->imageindex!=-1) {
printf("(%d,%d), ",current->domindex,current->imageindex);
current=current->next;
}
}
printf("]\n");
printf("edge (2,3) exists? -- %u\n", exists_edge(2,3,sam));
printf("edge (1,3) exists? -- %u\n", exists_edge(1,3,sam));
free(sam);
}*/
