#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "linear_algebra.h"
 
#define PI 4*atan(1)

void RK4(double CanonicalInitialCondition(double x), double x0, double xf, double dx, double dt, int J, int N){
	int i;
	struct vector *k1 = vector_create(J);
	struct vector *k2 = vector_create(J);
	struct vector *k3 = vector_create(J);
	struct vector *k4 = vector_create(J);

	struct vector *x = vector_create(J); //x-axis
	for (i = 0; i < x->n; i++){
		x->vals[i] = i * dx;
	}
	struct vector *ic = vector_create(J); //initial condition
	for (i = 0; i < ic->n; i++){
		ic->vals[i] = CanonicalInitialCondition(x->vals[i]);	
	}
	struct vector *f = vector_create(J); //dummy vector
	/////////////////////////////////////////////////////
	FILE *fout1 = NULL;   //x-axis for plotting
	fout1 = fopen("results_xaxis.txt", "w");
	for(i=0; i<x->n;i++){
		fprintf(fout1, "%g\n", x->vals[i]);
	}
	fclose(fout1);
	///////////////////////////////////////////////////////
	FILE *fout2 = NULL;  //print initial state: 
	fout2 = fopen("results_initial.txt", "w");
	for(i = 0; i < J; i++){
		fprintf(fout2, "%g\n", ic->vals[i]);
	}
	fclose(fout2);	
	///////////////////////////////////////////////////////
	//FILE *fout4 = NULL; //Interdiate States (in loop--we certainly do not want all of them)  
	//fout4 = fopen("results_intermediate.txt", "w");	
	////////////////////////////////////////////////////////

	//Start clock
	double time_naught = WTime();
	
	for(i = 0; i < N; i++){
		// k1 = f(u_i)
		rhs(ic,k1,dx);
		// k2 = f(U_i + dt/2*k1)
		vector_add(1.,ic,dt/2,k1,f);
		rhs(f,k2,dx);
		// k3 = f(U_i + dt/2*k2)
		vector_add(1.,ic,dt/2,k2,f);
		rhs(f,k3,dx);
		// k4 = f(U_i + dt*k3)
		vector_add(1.,ic,dt,k3,f);
		rhs(f,k4,dx);

		for(int j = 0; j < J; j++) {
			VEC(ic, j) = VEC(ic, j) + dt/6 * (VEC(k1,j)+ 2*VEC(k2,j)+ 2*VEC(k3, j)+ VEC(k4,j));
		}
		/*
		//Print (just for graphing--cutting this down will speed up code)
		if(i % 20 == 0){
		for(int k=0; k<J-1; k++){
			fprintf(fout4, "%.17le\t", ic->vals[k]);
		}
		fprintf(fout4,"\n");
		}
		*/
	}
	
	//fclose(fout4);

	double time_final = WTime();
	printf("Time = %g\n", time_final - time_naught);
	

	//prints final state: will duplicate intermediate data if user prints all
	FILE *fout3 = NULL;   
	fout3 = fopen("results_final.txt", "w");
	for(i = 0; i < J-1; i++){
		fprintf(fout3, "%g\n", ic->vals[i]);
	}
	fclose(fout3);




	vector_destroy(k1);
	k1 = NULL;
	vector_destroy(k2);
	k2 = NULL;
	vector_destroy(k3);
	k3 = NULL;
	vector_destroy(k4);
	k4 = NULL;
	vector_destroy(ic);
	ic = NULL;
	vector_destroy(f);
	f = NULL;
}




