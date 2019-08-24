#include <stdio.h>
#include <math.h>

int    n=3;
double f[]={6.0, 9.0, 2.0, 5.0};
double x[]={1.0, 2.0, 3.0, 4.0};
double c[4];

/******************************************************************
 * Procedure to compute Newton polynomial in forward difference.
 * Input:
 *  x[]: (n+1) parametric values,
 *  f[]: (n+1) function values,
 *  n: last index of samples, (n+1) samples
 * Output:
 *  c[]: coefficients of Newton's polynomial, forward divided difference
 *  
 */
void newton(double x[], double f[], double c[], int n)
{
  int    i, j, k;

  //Order 0 divided differences
  for(j=0;j<=n;j++)
    c[j] = f[j];
  //Order 1 to n divided differences
  for(k=1;k<=n;k++){
    for(i=n;i>=k;i--)
      c[i] = (c[i]-c[i-1])/(x[i]-x[i-k]);
  }
  for(k=0;k<=n;k++)
    fprintf(stderr, "C[%d]=%lf\n",k, c[k]);
}


/*--------------------------------------------------------
 * Procedure to evaluate Newton's polynomial.
 *  c[]: coefficients of Newton polynomial.
 *  x[]: parametric values,
 *  t: interpolation position inside range x[0]~x[n],
 *  n: degree of Newton's polynomial
 */
double Horner(double t, double c[], double x[], int n)
{
  double sum;
  int    i;

  sum = c[n];
  for(i=n-1;i>=0;i--)
    sum = sum*(t-x[i]) + c[i];
}


/*----------------------------------------------------------
 * The main procedure
 */
void main(int argc, char **argv)
{
  int     i;
  double  p[4];

  //Compute Newton's polynomial.
  newton(x, f, c, n);
  //Verify the interpolation conditions
  for(i=0;i<=n;i++){
    p[i] = Horner(x[i], c, x, n);
    fprintf(stderr,"i=%d, sample= %lf, interpolated= %lf\n", i, f[i], p[i]);
  }
}
