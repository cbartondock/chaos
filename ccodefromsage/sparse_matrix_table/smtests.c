#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "./smtable.c"

int main(char* args) {
    int grid=25;
    unsigned char m[grid][grid];
    for(int i=0; i < grid; i++) {
        for(int j=0; j<grid; j++) {
            m[i][j]=1;
        }
    }
    for(int i=0; i < grid; i++) {
        for(int j=0; j<grid; j++) {
            printf("%u  ",m[i][j]);
        }
        printf("\n");
    }
    sparse_adjacency_matrix* sm = initialize_sparse_matrix(grid, 2, 2,-2.5,-2.5,2.5,2.5,m);
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
