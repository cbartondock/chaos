#include <math.h>
#include <quadmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

__float128 bucketsum(__float128 *l, int num) {
    __float128 max = FLT128_MAX;
    __float128 min = FLT128_MIN;
    for(int i=0; i < num; i++) {
        if(fabsq(l[i]) > fabsq(max)) {
            max = l[i];
        }
        if(fabsq(l[i]) < fabsq(min)) {
            min = l[i];
        }
    }
    printf("max: %Lf, min: %Lf\n", (long double)max,(long double)min);
    int max_exp = (int)ceilq((log(fabsq(max))/logq(2.Q)));
    int min_exp = fmaxq(-100,(int)floorq((logq(fabsq(min))/logq(2.Q))));
    printf("max_exp: %d, min_exp: %d\n",max_exp,min_exp);
    int len = max_exp - min_exp;
    __float128 buckets[len]; 
    for(int m=0; m < len; m++) {
        buckets[m]=0.Q;
    }
    int exp;
    __float128 bucktop;
    __float128 buckbot;
    for(int i=0; i < num; i++) {
        for(int b=0; b < len; b++) {
            exp = min_exp + b;
            bucktop = powq(2.Q,exp+1.Q);
            buckbot = powq(2.Q,exp);
            if(fabsq(l[i]) <= bucktop && fabsq(l[i]) > buckbot) {
                buckets[b] += l[i];
                int k=b;
                if(l[i]>0.Q) {
                    while(fabs(buckets[k]) > powq(2.Q, min_exp+k+1.Q) && k < len-1.Q) {
                        buckets[k+1] += buckets[k];
                        buckets[k]=0.Q;
                        k++;
                    }
                }
                else if(l[i]<0.Q) {
                    while(fabsq(buckets[k]) < powq(2.Q, min_exp+k) && k > 0) {
                        buckets[k-1] += buckets[k];
                        buckets[k]=0.Q;
                        k--;
                    }
                }
                break;
            }
        }
    }
    __float128 result = 0.Q;
    for(int m=0; m < len; m++) {
        result+=buckets[m];
    }
    return result;
}

__float128 smod(__float128 x, __float128 y) {
    int div = (int)(x/y);
    __float128 result = x - div*y;
    result = (result < 0.Q) ? result + y : result;
    return result;
}

/*
int main() {
    long double *list = (long double *) malloc(sizeof(long double)*5);
    list[0] = 0.; list[1]=-.1; list[2] =.4; list[3] =-.401; list [4]= .33;
    printf("proper sum: %Lf\n",list[0]+list[1]+list[2]+list[3]+list[4]);
    printf("bucket sum: %Lf\n", bucketsum(list, 5));
}*/
