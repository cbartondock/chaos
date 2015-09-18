#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../usefulfunctions/functions.c"
long double weight(long double t) {
    if(t<=0 || t>=1) { return 0.;}
    return exp(1/(t*(t-1)));
}

int main() {
    FILE *out = fopen("outputs/sumdiffs.txt", "w");
    for(int i=3; i<=8; i++) {
        int time = pow(10,i);
        long double wsuml = 0.;
        double wsum = 0.;
        double bsum;
        long double * ws = (long double *) malloc(sizeof(long double)*time);
        for(int t=0; t < time; t++) {
            ws[t] = weight((long double)t/(long double) time);
            wsum += weight((long double)t/(long double) time);
            wsuml += weight((long double)t/(long double)time);
        }
        bsum = bucketsum(ws,time);
        fprintf(out, "time is: %u\n", time); 
        fprintf(out, "wsum normal: %.50f\n", wsum);
        fprintf(out, "wsum long: %.50Lf\n", wsuml);
        fprintf(out, "wsum buckets: %.50f\n", bsum);
        fprintf(out, "digit of diff bd: %Lf\n", (long double)-log(fabs(bsum-wsum))/log(10));
        fprintf(out, "digit of diff bl: %Lf\n", (long double)-log(fabs(bsum-wsuml))/log(10));
        fprintf(out, "digit of diff ld: %Lf\n", (long double)-log(fabs(wsum-wsuml))/log(10));
        fprintf(out, "%% error bl: %Lf\n", (long double) -log(fabs(bsum - wsuml)/fabs(wsuml))/log(10));
        fprintf(out, "%% error ld: %Lf\n", (long double) -log(fabs(wsuml-wsum)/fabs(wsuml))/log(10));
    }
}
