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
    double leastx;
    double leasty;
    double deltax;
    double deltay;
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
sparse_adjacency_matrix* initialize_sparse_matrix(int grid,int chmap,int top,double leastx, double leasty, double deltax, double deltay, unsigned char (*m)[grid]) {
    int domnumber=0;
    for(int i=0; i<grid; i++) {
        for(int j=0; j<grid; j++) {
            domnumber+=m[i][j];
        }
    }
    sparse_adjacency_matrix* adjmatrix;
    adjmatrix = (sparse_adjacency_matrix*) malloc(sizeof(sparse_adjacency_matrix));
    adjmatrix->adjacency_lists = (adj_element**) malloc(sizeof(adj_element*)*domnumber);
    adjmatrix->grid = grid;
    adjmatrix->leastx = leastx;
    adjmatrix->leasty = leasty;
    adjmatrix->deltax = deltay;
    adjmatrix->deltay = deltay;
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
    int i,j, domindex, imageindex;
    double corners[4][2];
    double newcorners[4][2];
    double edges[4][2][2];
    double edge[2][2];
    double diff[2];
    double tiv,tjv, ceptx, cepty;
    int iceptx, icepty,iceptx2,icepty2;
    for(int b=0; b<domnumber; b++) {
        domindex = (adjmatrix->adjacency_lists)[b]->domindex;
        i=domindex/grid;
        j=domindex%grid;

        corners[0][0] = leastx + j*deltax; corners[0][1] = leasty + i*deltay; //tl
        corners[1][0] = leastx + (j+1)* deltax; corners[1][1] = leasty + i*deltay; //tr
        corners[2][0] = leastx + j*deltax; corners[2][1] = leasty + (i+1)*deltay; //bl
        corners[3][0] = leastx + (j+1)*deltax; corners[3][1] = leasty + (i+1)*deltay; //br
        for(int c=0; c<4; c++) {
            //double x = corners[c][0];
            //double y = corners[c][1];
            /*if(chmap==1){
              corners[c][1] = smod(sin(x+y)+y,2*M_PI);
              corners[c][0] = smod(x+y,2*M_PI);
              }*/
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

            int ivmin = (int)ceil((min(edge[0][1],edge[1][1])-leasty)/deltay);
            int ivmax = (int)floor((max(edge[0][1],edge[1][1])-leasty)/deltay);
            int jvmin = (int)ceil((min(edge[0][0],edge[1][0])-leastx)/deltax);
            int jvmax = (int)floor((max(edge[0][0],edge[1][0])-leastx)/deltax);

            for(int iv=ivmin; iv<=ivmax; iv++) {
                tiv = (leasty+iv*deltay-edge[0][1])/(diff[1]);
                ceptx = edge[0][0] + diff[0]*tiv;
                cepty = edge[0][1] + diff[1]*tiv;
                iceptx = (int)round((ceptx-leastx)/deltax);
                icepty = (int)round((cepty-leasty)/deltay);
                iceptx2= iceptx-1;
                icepty2= icepty-1;
                if(top==2){
                    //cylinder
                    iceptx = smod(iceptx,grid);
                    iceptx2 = smod(iceptx2,grid);
                }
                imageindex = iceptx + grid*icepty;
                if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                    add_edge(domindex,imageindex,adjmatrix);
                }
                imageindex = iceptx + grid*icepty2;
                if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                    add_edge(domindex,imageindex,adjmatrix);
                }
            }
            for(int jv=jvmin; jv<=jvmax; jv++) {
                tjv = (leastx+jv*deltax-edge[0][0])/(diff[0]);
                ceptx = edge[0][0] + diff[0]*tjv;
                cepty = edge[0][1] + diff[1]*tjv;
                iceptx = (int)round((ceptx-leastx)/deltax);
                icepty = (int)round((cepty-leasty)/deltay);
                iceptx2= iceptx-1;
                icepty2= icepty-1;
                if(top==2){
                    //cylinder
                    iceptx = smod(iceptx,grid);
                    iceptx2 = smod(iceptx2,grid);
                }
                imageindex = iceptx + grid*icepty;
                if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                    add_edge(domindex,imageindex,adjmatrix);
                }
                imageindex = iceptx2 + grid*icepty;
                if(exists_vertex(imageindex,adjmatrix) && !exists_edge(domindex,imageindex,adjmatrix)) {
                    add_edge(domindex,imageindex,adjmatrix);
                }

            }
            if(jvmax < jvmin && ivmax < ivmin) {
                m[i][j] = (1<<4)|m[i][j];
            }
        }
    }
    return adjmatrix;
}

void free_matrix(sparse_adjacency_matrix* matrix) {
    adj_element* current;
    adj_element* next;
    for(int i=0; i<matrix->grid; i++) {
        current = (matrix->adjacency_lists)[i]->next;
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
int grid(sparse_adjacency_matrix* matrix) {
    return matrix->grid;
}
int domnumber(sparse_adjacency_matrix* matrix) {
    return matrix->domnumber;
}

/*sparse_adjacency_matrix* quadruple_precision(sparse_adjacency_matrix* oldmatrix) {
  sparse_adjacency_matrix* newmatrix;
  newmatrix = (sparse_adjacency_matrix*) malloc(sizeof(sparse_adjacency_matrix));
  newmatrix->adjacency_lists = (adj_element**) malloc(sizeof(adj_element*)*(oldmatrix->domnumber)*4);
  newmatrix->grid = 2*(oldmatrix->grid);
  newmatrix->domnumber = 4*(oldmatrix->domnumber);
  newmatrix->leastx = oldmatrix->leastx;
  newmatrix->leasty = oldmatrix->leasty;
  newmatrix->deltax = oldmatrix->deltax/2.;
  newmatrix->deltay = oldmatrix->deltay/2.;
  adj_element* newelement;
  int oldgrid = oldmatrix->grid;
  int newgrid = newmatrix->grid;
  int olddomindex;
  int domindices[4]={0,0,0,0};
  for(int n=0; n < oldmatrix->domnumber; n++) {
  olddomindex = (oldmatrix->adjacency_lists)[n]->domindex;
  domindices[0] = 2*(olddomindex % oldgrid) + 2*(olddomindex/oldgrid)*(newgrid);
  domindices[1] = 2*(olddomindex % oldgrid) + 1 + 2*(olddomindex/oldgrid)*(newgrid);
  domindices[2] = 2*(olddomindex % oldgrid) + (2*(olddomindex/oldgrid) + 1)*(newgrid);
  domindices[3] = 2*(olddomindex % oldgrid) + 1 + (2*(olddomindex/oldgrid) + 1)*(newgrid);
  for(int k=0; k<4; k++) {
  newelement = (adj_element*) malloc(sizeof(adj_element));
  newelement->domindex = domindices[k];
  newelement->imageindex = -1;
  newelement->imagenumber = 0;
  newelement->next = newelement;
  newelement->prev = newelement;
  (newmatrix->adjacency_lists)[4*n+k]=newelement;
  }
  }
  int i,j,newdomindex,newimageindex;
  double x, y, xn, yn;
  adj_element* current;
  for(int n=0; n < oldmatrix->domnumber; n++) {
  olddomindex = (oldmatrix->adjacency_lists)[n]->domindex;
  for(int m=4*n; m < 4*n+4; m++) {
  newdomindex = ((newmatrix->adjacency_lists)[m]->domindex);
  i = newdomindex/(newmatrix->grid);
  j = newdomindex%(newmatrix->grid);
  for(int s=0; s <=newmatrix->numsamples; s++){
  for(int s2=0; s2 <=newmatrix->numsamples; s2++) {
  x = newmatrix->leastx + j*(newmatrix->deltax) + newmatrix->deltax*((double)s/(double)newmatrix->numsamples);
  y = newmatrix->leasty + i*(newmatrix->deltay) + newmatrix->deltay*((double)s2/(double)newmatrix->numsamples);
  xn = (*(newmatrix->fx))(x,y);
  yn = (*(newmatrix->fy))(x,y);
  newimageindex = (int)floor((xn-(newmatrix->leastx))/(newmatrix->deltax)) + (newmatrix->grid)*(int)floor((yn-(newmatrix->leasty))/newmatrix->deltay);
  if(exists_vertex(newimageindex,newmatrix) && !exists_edge(newdomindex,newimageindex,newmatrix)) {
  add_edge(newdomindex,newimageindex,newmatrix);
  }
  }
  }
  }
  }
  return newmatrix;
  }*/

double henonx(double x, double y) {
    return 1.4 - x*x + .3*y;
}
double henony(double x, double y) {
    return x;
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
