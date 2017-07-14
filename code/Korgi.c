#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "linear_algebra.h"
 
#define PI 4*atan(1)

double CanonicalInitialCondition(double x){
    double velocity1 = 128.;
    double amplitude1 = velocity1 / 2.;
    double wave_number1 = sqrt(velocity1) / 2.;
    double velocity2 = 64.;
    double amplitude2 = velocity2 / 2.;
    double wave_number2 = sqrt(velocity2) / 2.;    
    return (amplitude1 * (1. / cosh(wave_number1 * x - 10.)) * (1. / cosh(wave_number1 * x - 10.)))  + (amplitude2 * (1. / cosh(wave_number2 * x - 20.)) * (1. / cosh(wave_number2 * x - 20.)));
}

  

int main(){
	double x0 = 0.;						//initial x position
	double xf = 50.;					//final x position
	double dt = 0.00001; 				//temporal step size
	int J = 1000;						//space steps
	int N = 10000;						//time steps
	double dx = (xf - x0) / (J - 1.);	//spatial step size
	RK4(CanonicalInitialCondition, x0, xf, dx, dt, J, N);
}