#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <quadmath.h>
/*
 * Performs a 4th order 2-dimensional Runge Kutta, 
 * with the initial conditions stored in p.
 * Doesn't return a value, it modifies p.
 * Runs from 0 to tmax with timestep h.
 */
__float128 ddpx(__float128 x, __float128 y, __float128 t) {
    return y;
}
__float128 ddpy(__float128 x, __float128 y, __float128 t) {
    return 2.5Q*cos(t)-sin(x)-0.2Q*y;
}



void rk4(__float128 p[2], __float128 tmax, __float128 h) {
    __float128 t = 0;
    __float128 (*f)(__float128,__float128,__float128) = &ddpx;
    __float128 (*g)(__float128,__float128,__float128) = &ddpy;
    __float128 x = p[0];
    __float128 y = p[1];
    __float128 k1,k2,k3,k4,l1,l2,l3,l4;
    while(t <= tmax) {
        k1 = (*f)(x,y,t);
        l1 = (*g)(x,y,t);
        k2 = (*f)(x + 0.5Q*h*k1, y + 0.5Q*h*l1, t + 0.5Q*h);
        l2 = (*g)(x + 0.5Q*h*k1, y + 0.5Q*h*l1, t + 0.5Q*h);
        k3 = (*f)(x + 0.5Q*h*k2, y + 0.5Q*h*l2, t + 0.5Q*h);
        l3 = (*g)(x + 0.5Q*h*k2, y + 0.5Q*h*l2, t + 0.5Q*h);
        k4 = (*f)(x + h*k3, y + h*l3, t + h);
        l4 = (*g)(x + h*k3, y + h*l3, t + h);
        x += (h/6.Q)*(k1 + 2.Q*k2 + 2.Q*k3 + k4);
        y += (h/6.Q)*(l1 + 2.Q*l2 + 2.Q*l3 + l4);
        t += h;
    }
    p[0]=x;
    p[1]=y;
}
/*
int main() {
    printf("running main\n");
    __float128 initial[2];
    initial[0] = .38;
    initial[1] = 5Q; 
    rk4(initial,5Q00*M_PI,.1);
    printf("theta: %f, Dtheta: %f\n",initial[0],initial[1]);
}
*/

