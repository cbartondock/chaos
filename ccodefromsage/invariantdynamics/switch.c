#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../usefulfunctions/functions.c"
void switchcount() {
    double x = .01;
    double y = .01;
    double xn, yn;
    double average_x=0;
    double average_y=0;
    double last_average_x=0;
    double last_average_y=0;
    double diff;
    int totaltime = 1000000000;
    int averagetime = 10000;
    double thresh=.3;
    int switchcount = 0;
    for(int t=1; t<= totaltime; t++) {
        average_x+=x;
        average_y+=y;
        if(t%averagetime == 0) {
            average_x/=averagetime; average_y/=averagetime;
            diff = (last_average_x-average_x)*(last_average_x-average_x) + (last_average_y-average_y)*(last_average_y-average_y);
            if(t> averagetime && diff > thresh) {

                printf("switchcount: %u\n",switchcount);
                printf("last average is (%f, %f)\n",last_average_x, last_average_y);
                printf("average is (%f, %f)\n\n",average_x,average_y);
                switchcount+=1;
            }

            last_average_x=average_x;
            last_average_y=average_y;
            average_x=0;
            average_y=0;
        }

        xn = smod(x+y, 2*M_PI);
        yn = smod(1.4*sin(x+y)+y, 2*M_PI);
        x=xn;
        y=yn;
    }
    printf("========================================\n\n");
    printf("switch probability is: %.20f\n",(double)switchcount/(double)totaltime);

}




int main() {
    switchcount();
    return 0;
}
