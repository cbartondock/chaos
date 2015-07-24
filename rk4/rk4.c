#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
 * Performs a 4th order 2-dimensional Runge Kutta, 
 * with the initial conditions stored in p.
 * Doesn't return a value, it modifies p.
 * Runs from 0 to tmax with timestep h.
 */
double ddpx(double x, double y, double t) {
    return y;
}
double ddpy(double x, double y, double t) {
    return 1.8*cos(t)-sin(x)-.2*y;
}



void rk4(double p[2], double tmax, double h) {
    double t = 0;
    double (*f)(double,double,double) = &ddpx;
    double (*g)(double,double,double) = &ddpy;
    double x = p[0];
    double y = p[1];
    double k1,k2,k3,k4,l1,l2,l3,l4;
    while(t <= tmax) {
        k1 = (*f)(x,y,t);
        l1 = (*g)(x,y,t);
        k2 = (*f)(x + 0.5*h*k1, y + 0.5*h*l1, t + 0.5*h);
        l2 = (*g)(x + 0.5*h*k1, y + 0.5*h*l1, t + 0.5*h);
        k3 = (*f)(x + 0.5*h*k2, y + 0.5*h*l2, t + 0.5*h);
        l3 = (*g)(x + 0.5*h*k2, y + 0.5*h*l2, t + 0.5*h);
        k4 = (*f)(x + h*k3, y + h*l3, t + h);
        l4 = (*g)(x + h*k3, y + h*l3, t + h);
        x += (h/6.)*(k1 + 2*k2 + 2*k3 + k4);
        y += (h/6.)*(l1 + 2*l2 + 2*l3 + l4);
        t += h;
    }
    p[0]=x;
    p[1]=y;
}
/*
int main() {
    printf("running main\n");
    double initial[2];
    initial[0] = .38;
    initial[1] = 5; 
    rk4(initial,500*M_PI,.1);
    printf("theta: %f, Dtheta: %f\n",initial[0],initial[1]);
}
*/
