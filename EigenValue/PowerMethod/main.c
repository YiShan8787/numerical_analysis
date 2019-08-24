/*************************************************************
* In this sample program. the power method and inverse power
* method are used to compute the max and min eigenvalues of
* a matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "definition.h"
#include "vec_mtx.h"

/*--- Declare the target matrix. ---*/
double   **A;


/*--- Declare the max and min eigenvalues and their associated
      eigenvectors.
	  */
double  minEigenVal, maxEigenVal, *minEigenVec, *maxEigenVec;




/*----------------------------------------------------
 * The main procedure
 */
int main(int argc, char **argv)
{
	int  i, j, n;
	int  k;
	double   t;

	n = 3;

    A = alloc_mtx(n);

	/*--- Initialize A[]. ---*/
    for(i=0;i<n;i++)
		for(j=i;j<n;j++) {
//			k = rand() % 100;
			A[i][j] = A[j][i] = (double) i+j+1.0;
 //   		A[i][j] = A[j][i] = 1.0/(i+j+1.0);
		}

	// Compute the dominant eigenvalue and its
	// associated eigenvector.
	maxEigenVec = alloc_vec(n);
    minEigenVec = alloc_vec(n);
	power_method(A, maxEigenVec, &maxEigenVal, n);
    inverse_power_method(A, minEigenVec, &minEigenVal, n);

	fprintf(stderr,"A[][]=\n");
    for(i=0;i<n;i++){
		fprintf(stderr,"\n");
		for(j=0;j<n;j++)
			fprintf(stderr,"%lf ", A[i][j]);
	}
	fprintf(stderr,"\n--------------------------\n");

	fprintf(stderr,"Max eigenvalue= %lf\n", maxEigenVal);
	fprintf(stderr,"The max eigen vector=\n");
    for(i=0;i<n;i++)
      fprintf(stderr," %lf ", maxEigenVec[i]);
	fprintf(stderr,"\n");
	// print out min eigen value & vector
	fprintf(stderr,"Min eigenvalue= %lf\n", minEigenVal);
	fprintf(stderr,"The min eigen vector=\n");
    for(i=0;i<n;i++)
      fprintf(stderr," %lf ", minEigenVec[i]);
	fprintf(stderr,"\n");
    t = inner_product(minEigenVec, maxEigenVec, n);
	fprintf(stderr," <Vmax, Vmin>= %lf\n", t);
	fprintf(stderr," Condition number = %lf\n", fabs(maxEigenVal/minEigenVal));
	getchar();
}
