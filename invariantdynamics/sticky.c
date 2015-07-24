#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../usefulfunctions/functions.c"

int sum_matrix(int grid, int (*matrix)[grid]) {
    int result=0;
    for(int i=0; i < grid; i++) {
        for(int j=0; j < grid; j++) {
            result+=matrix[i][j];
        }
    }
    return result;
}

void rp_window(int windowwidth, double ix, double iy, double epsilon, int (*m)[windowwidth], double xlist[windowwidth], double ylist[windowwidth]) {
    double x=ix;
    double y=iy; 
    double xn,yn,xn2,yn2;
    for(int i=0; i < windowwidth; i++) {
        for(int j=0; j < windowwidth; j++) {
            m[i][j]=0;
        }
    }
    for(int t=0; t< windowwidth; t++) {
        double x2=ix;
        double y2=iy;
        for(int t2=0; t2<windowwidth; t2++) {
            if(fabs(x2-x)<epsilon && fabs(y2-y) <epsilon) {
                m[windowwidth-1-t2%windowwidth][t%windowwidth]+=1;
            }
            xn2 = smod(x2+y2, 2*M_PI);
            yn2 = smod(1.4*sin(x2+y2)+y2, 2*M_PI);
            x2=xn2;
            y2=yn2;
        }
        xlist[t] = x;
        ylist[t] = y;
        xn = smod(x+y, 2*M_PI);
        yn = smod(1.4*sin(x+y)+y, 2*M_PI);
        x=xn;
        y=yn;
    }

    for(int i=0; i < 10; i++) {
        for(int j=0; j < 10; j++) {
            printf("%u ",m[i][j]);
        }
        printf("\n");
    }
}

void stickiness(int windows,double ix, double iy, double epsilon, double rrs[windows]) {

    FILE *points_out = fopen("outputs/stickypoints.txt","w");
    if (points_out == NULL) {
        printf("file error");
        exit(1);
    }
    printf("\n");

    double x = ix;
    double y = iy;
    double xn, yn;
    int windowwidth = 5000;
    int total = windowwidth*windowwidth;
    int windowlap = 4500;
    int (*m)[windowwidth]= (int(*)[windowwidth]) malloc(sizeof(int[windowwidth])*windowwidth);
    for(int i=0; i < windowwidth; i++) {
        for(int j=0; j < windowwidth; j++) {
            m[i][j]=0;
        }
    }
    int in_sticky = 0;
    double xlist[windowwidth];
    double ylist[windowwidth];
    for(int w=0; w<windows; w++) {
        rp_window(windowwidth,x,y,epsilon,m, xlist, ylist);
        double rr = (double)sum_matrix(windowwidth,m)/(double)total;
        printf("w: %u, rr: %f\n",w,rr);
        rrs[w] = rr;
        if(rr>0.05) {
            printf("xlist[0]: %f\n", xlist[0]);

            if(!in_sticky) {
                printf("[");
                for(int i=0; i<windowwidth; i++) {
                    printf("(%f, %f), ",xlist[i],ylist[i]);
                    fprintf(points_out, "(%f, %f), ",xlist[i],ylist[i]);
                }
                printf("]\n");
                in_sticky=1;
            }
            if(in_sticky) {
                printf("[");
                for(int i=windowlap; i < windowwidth; i++) {
                    printf("(%f, %f), ",xlist[i],ylist[i]);
                    fprintf(points_out, "(%f, %f), ",xlist[i],ylist[i]);
                }
                printf("]\n");
            }
        } else {
            in_sticky = 0;
        }

        for(int t=0; t<windowwidth-windowlap; t++) {
            xn = smod(x+y,2*M_PI);
            yn = smod(1.4*sin(x+y)+y,2*M_PI);
            x=xn;
            y=yn;
        }
    }
    fclose(points_out);
}


int main() {
    int windowwidth = 10; 
    double xlist[windowwidth];
    double ylist[windowwidth];
    int (*m)[windowwidth] = (int(*)[windowwidth]) malloc(sizeof(int[windowwidth])*windowwidth);
    for(int s=0; s< windowwidth; s++) {
        for(int s2 =0; s2<windowwidth;s2++) {
            m[s][s2]=0;
        }
    }
    rp_window(windowwidth,M_PI,1.2,.5,m,xlist,ylist);
    /*int windows = 3;
      double rrs[windows];
      stickiness(windows,5.14,2.14,.5,rrs);
      printf("rr's: [");
      for(int w=0; w<windows; w++) {
      printf("%f, ",rrs[w]);
      }
      printf("]\n");*/
}
