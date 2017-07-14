
#include "linear_algebra.h"

#include <omp.h>

// ----------------------------------------------------------------------
// vector_add
//
// calculate z = x + y, for the vectors x, y, z
void vector_add(double a, struct vector *x, double b, struct vector *y, struct vector *z){
	// can only add vectors of equal length
	assert(x->n == y->n);
	assert(x->n == z->n);
	#pragma omp for
	for (int i = 0; i < x->n; i++) {
		VEC(z, i) = a * VEC(x, i) + b * VEC(y, i);
	}
}
