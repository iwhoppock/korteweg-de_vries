#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "linear_algebra.h"
 
#define PI 4*atan(1)


void rhs(struct vector *x, struct vector *v, double dx){
	int n = v->n;
	//Initial Step
	int	i = 0;
	VEC(v, i) = -6 * VEC(x,i) * (VEC(x,i+1) - VEC(x,n-1)) / 2 / dx;
	VEC(v, i) += - (VEC(x,i+2) - 2 * VEC(x, i+1) + 2 * VEC(x,i+n-1) - VEC(x,i-2+n)) / 2 / (dx * dx * dx);
	//First Step
	i = 1;
	VEC(v, i) = -6 * VEC(x,i) * (VEC(x,i+1) - VEC(x,i-1)) / 2 / dx;
	VEC(v, i) += - (VEC(x,i+2) - 2 * VEC(x, i+1) + 2 * VEC(x,i-1) - VEC(x,i-2+n)) / 2 / (dx * dx * dx);
	//Boundary Conditions
	#pragma omp for
	for (i = 2; i < n - 2; i++){
		VEC(v, i) = -6 * VEC(x,i) * (VEC(x,i+1) - VEC(x,i-1)) / 2 / dx;
		VEC(v, i) += -(VEC(x,i+2) - 2 * VEC(x, i+1) + 2 * VEC(x,i-1) - VEC(x,i-2)) / 2 / (dx * dx * dx);
	}
	//Outside of Boundary Point One
  	i = n-2;
  	VEC(v, i) = -6 * VEC(x,i) * (VEC(x,i+1) - VEC(x,i-1)) / 2 / dx;
  	VEC(v, i) += - (VEC(x,i+2-n) - 2 * VEC(x, i+1) + 2 * VEC(x,i-1) - VEC(x,i-2)) / 2 / (dx * dx * dx);
	//Outside of Boundary Point Two
  	i = n-1;
  	VEC(v, i) = -6 * VEC(x,i) * (VEC(x,i+1-n) - VEC(x,i-1)) / 2 / dx;
  	VEC(v, i) += - (VEC(x,i+2-n) - 2 * VEC(x, i+1-n) + 2 * VEC(x,i-1) - VEC(x,i-2)) / 2 / (dx * dx * dx);
}





