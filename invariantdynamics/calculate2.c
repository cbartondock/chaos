#define __USE_GUU
#include <stdio.h>
#include <math.h>
#define pi M_PIl

static long double twopi = 6.283185307179586476925286;

long double mymods(long double x);
long double calculate_x(long double x, long double y);
long double calculate_y(long double x, long double y);
int process(long double x, long double y);
long double w(long double t);
int count();

int main(){
    long double step;
    int i;
    int j;
    long double inix;
    long double iniy;
    int times;
    int result;

    int count[22];

    for(i = 0; i <22; i++){
        count[i]=0;
    }

    inix = 0;
    iniy = 0;

    times = 20;

    step = (twopi - inix) / (long double) times;

    for(i = 0; i < times; i++){
        for(j = 0; j < times; j++){
            result = process(inix, iniy);

            if(result < 0){
                result = 0;
            }

            if(result > 20){
                result = 20;
            }

            count[result] = count[result] + 1;
            inix = inix + step;
        }

        inix = 0;
        iniy = iniy + step;
    }
    for(i = 0; i < 21; i ++){
        printf("the count of %d  is %d \n" ,i , count[i]);

    }



}


int process(long double x, long double y){
    int i;	
    /* times is how many iterations we want */
    int const times = 1000;
    long double x_cm;
    long double y_cm;
    long double x_cm2;
    long double y_cm2;
    long double temx;
    long double temy;
    long double temx2;
    long double temy2;
    long double t;
    long double weight_sum;
    long double diff_x;
    long double diff_y;
    long double norm;
    int result;

    x_cm = 0;
    y_cm = 0;
    x_cm2 = 0;
    y_cm2 = 0;
    weight_sum = 0;

    printf("Initially: x is %Lf, y is %Lf \n", x, y);
    //printf("%Lf %Lf ", x, y);

    for(i = 0; i < times; i++){
        t = (long double) i / (long double) times;

        x_cm = x_cm + sin(x) * w(t);
        y_cm = y_cm + sin(y) * w(t);

        weight_sum = weight_sum + w(t);

        x = calculate_x(x, y);
        y = calculate_y(x, y);

        //printf("x is %Lf, y is %Lf \n", x ,y);
    }

    //printf("x now is %Lf, y is %Lf \n", x, y);
    temx = x_cm / weight_sum;
    temy = y_cm / weight_sum;
    //printf("x cumulative is %Lf, y is %Lf \n", x_cm, y_cm);
    //printf("x average is %.20Lf, y is %.20Lf \n", temx, temy);


    weight_sum = 0;

    for(i = 0; i < times; i++){
        t = (long double) i / (long double) times;
        x_cm2 = x_cm2 + sin(x) * w(t);
        y_cm2 = y_cm2 + sin(y) * w(t);

        weight_sum = weight_sum + w(t);

        x = calculate_x(x, y);
        y = calculate_y(x, y);

        //printf("x is %Lf, y is %Lf \n", x, y);
    }

    //printf("x now is %Lf, y now is %Lf \n",x ,y);
    temx2 = x_cm2 / weight_sum;
    temy2 = y_cm2 / weight_sum;
    //printf("x cumulative2 is %Lf, y is %Lf \n", x_cm2, y_cm2);
    //printf("x average is %.20Lf, y is %.20Lf \n", temx2, temy2);

    diff_x = temx - temx2;
    diff_y = temy - temy2;

    printf("difference of x average is %.20Lf\n", temx - temx2);
    printf("difference of y average is %.20Lf\n", temy - temy2);

    norm = sqrt(diff_x * diff_x + diff_y * diff_y );
    norm = -log(norm) / (long double) log(10);

    printf("norm is %.1Lf\n", norm);

    result = (int) norm;
    //printf("%d\n",result);

    printf("\n\n\n");

    if(result < 0){
        result = 0;
    }

    //printf("%d ", result);
    //printf("%f ", log(abs(diff_x)));
    //printf("%f \n", log(abs(diff_y)));
    return result;
}


long double calculate_x(long double x, long double y){
    long double x1;

    x1 = x + sin(y);
    x1 = mymods(x1);
    return x1;
}

long double calculate_y(long double x, long double y){
    long double y1;
    y1 = x + y;
    y1 = mymods(y1);
    return y1;
}

long double mymods ( long double a){

    while( a > twopi){
        a = a - twopi;
    }
    while(a < 0){
        a = a + twopi;
    }

    return a;
}

long double w (long double t){
    long double result;
    if((t == 0)||(t > 1)){
        return 0;
    }

    result = t * (t- 1);
    result = 1 / result;
    result = exp(result);
    return result; 
}

