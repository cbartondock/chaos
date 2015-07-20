#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "./smtable.c"

int main(char* args) {
    int grid=50;
    unsigned char m[grid][grid];
    for(int i=0; i < grid; i++) {
        for(int j=0; j<grid; j++) {
            m[i][j]=0;
        }
    } 
    m[-1+grid/2][-1+grid/2]=1;
    m[-1+grid/2][grid/2]=1;
    m[grid/2][-1+grid/2]=1;
    m[grid/2][grid/2]=1;
    for(int i=0; i < grid; i++) {
        for(int j=0; j<grid; j++) {
            printf("%u ", m[i][j]);
        }
        printf("\n");
    }

    sparse_adjacency_matrix* sm = initialize_sparse_matrix(grid,4,1,-3,-3,6./grid,6./grid,m);
    adj_element* current;
    printf("sam's edges are: [");
    for(int i=0; i <sm->domnumber; i++) {
        current = (sm->adjacency_lists)[i]->next;
        while(current->imageindex!=-1) {
            printf("(%d,%d), ",current->domindex,current->imageindex);
            current=current->next;
        }
    }
    printf("]\n");
    return 0;
}
