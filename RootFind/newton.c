/***********************************************************
 * This program demonstrate an example of root finding by using
 * Newton's method.
 */
#include <stdio.h>
#include <math.h>

#define EPSILON   0.000001

double   a, b;


/*---------------------------------------------------------------------------
 * The function fo the root-finding problem.
 */
double f(double x)
{
    return (x*x-5.0*x - 6.0);
}


/*--------------------------------------------------------------------------
 * derivative of the target function.
 */
double fx(double x)
{
  return(2.0*x-5.0);
}

/*--------------------------------------------------------------------------
 * The bisection method.
 *     a: the initial value of the root.
 */
double  newton(double a)
{
  int     i=1;
  double  err, xnew, xold;

  
  //Set up initial conditions.
  xold = a;
  xnew = xold - f(xold)/fx(xold);
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
