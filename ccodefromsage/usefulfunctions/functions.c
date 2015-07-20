#include <math.h>
#include <stdlib.h>
#include <stdio.h>
double smod(double x, double y) {
    int div = (int)(x/y);
    double result = x - div*y;
    result = (result < 0) ? result + y : result;
    return result;
}
/*int main() {
    printf("smod(.1,1) is: %f\n",smod(.1,1));
    printf("smod(-.1,1) is: %f\n",smod(-.1,1));
    printf("smod(1.1,1) is: %f\n",smod(1.1,1));
    printf("smod(-1.1,1) is: %f\n",smod(-1.1,1));
    printf("ceil(.1) is: %f\n",ceil(.1));
    printf("floor(.1) is: %f\n",floor(.1));
    printf("round(.1) is: %f\n",round(.1));
    printf("round(.6) is: %f\n",round(.6));

    printf("ceil(-.1) is: %f\n",ceil(-.1));
    printf("floor(-.1) is: %f\n",floor(-.1));
    printf("round(-.1) is: %f\n",round(-.1));
    printf("round(-.6) is: %f\n",round(-.6));
}*/
