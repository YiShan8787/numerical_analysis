/***********************************************************
 * This program demonstrate an example of root finding by using
 * Newton's method.
 */
#include <stdio.h>
#include <math.h>

#define EPSILON   0.000001

double   a, b;


/*---------------------------------------------------------------------------
 * The function f(x, y) for the root-finding problem.
 */
double f(double x, double y)
{
    return (x*x-5.0*y - 6.0);
}

/*---------------------------------------------------------------------------
 * The function g(x, y) for the root-finding problem.
 */
double g(double x, double y)
{
    return (3.0*x-y*y - 8.0);
}

/*-------------------------------------------------
 * Derivatives of f() and g().
 */
double fx(double x, double y)
{
	return 2.0*x;
}

double fy(double x, double y)
{
	return (-5.0);
}

double gx(double x, double y)
{
	return (3.0);
}

double gy(double x, double y)
{
	return (-2.0*y);
}




/*--------------------------------------------------------------------------
 * The bisection method.
 *     x, y: the initial value of the root.
 */
double  newton(double *x, double *y)
{
  int     i=1;
  double  err, xnew, ynew, xold, yold;
  double  h, k;
  double  fdx, fdy, gdx, gdy, Delta;
  double  fold, gold;

  
  //Set up initial conditions.
  xold = *x;
  yold = *y;

  fold = f(xold, yold);
  gold = g(xold, yold);
  fdx = fx(xold, yold);
  fdy = fy(xold, yold);
  gdx = gx(xold, yold);
  gdy = gy(xold, yold);

  Delta = fdx*gdy - fdy*gdx;

  xnew = xold - (;
  err = fabs(xnew - xold);

  fprintf(stderr, "     i            xn               error\n");
  fprintf(stderr,"------------------------------------------------------------\n");

  while(err>EPSILON){
   fprintf(stderr,"     %3d \t%lf \t%lf\n", i, xnew, err);

   xold = xnew; // Save current value.
   xnew = xold - f(xold)/fx(xold); // Compute new values.
   err = fabs(xnew - xold); // Compute the difference.
   i ++;
  }   
  fprintf(stderr,"     %3d \t%lf \t%lf\n", i, xnew, err);

  fprintf(stderr,"------------------------------------------------------------\n"); 
  return(xnew);
}


/*----------------------------------------------------------------------------
  * The main procedure.
  */
int main(int argc, char **argv)
{
   double    r;


   fprintf(stderr," Input the initial value x0= ");
   fscanf(stdin,"%lf", &r);
   r = newton(r);
   fprintf(stderr,"Root= %lf, error=%lf\n", r, f(r));
   fgetchar();
}
