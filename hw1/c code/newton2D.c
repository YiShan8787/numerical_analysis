/***********************************************************
 * This program demonstrate an example of root finding by using
 * Newton's method.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define EPSILON   0.000001

double   a, b;


/*---------------------------------------------------------------------------
 * The function f(x, y) for the root-finding problem.
 */
double f(double x, double y)
{
    return ((x-3)*(x-3)/25)+((y-4)*(y-4)/16)-1;//x*x-5.0*y - 6.0; ((x-3)*(x-3)/25)+((y-4)*(y-4)/16);
}

/*---------------------------------------------------------------------------
 * The function g(x, y) for the root-finding problem.
 */
double g(double x, double y)
{
    return ((x-3)*(x-3)/4)-((y-4)*(y-4)/9)-1;//3.0*x-y*y - 8.0; ((x-3)*(x-3)/4)+((y-4)*(y-4)/9)
}

/*-------------------------------------------------
 * Derivatives of f() and g().
 */
double fx(double x, double y)
{
	return (2.0*x-6)/25;//2.0*x; (2.0*x-6)/25
}

double fy(double x, double y)
{
	return (1.0*y-4)/8;//-5.0; (1.0*y-4)/8
}

double gx(double x, double y)
{
	return (1.0*x-3)/2;//3.0; (1.0*x-3)/2
}

double gy(double x, double y)
{
	return -(2.0*y-8)/9;//-2.0*y; (2.0*y-8)/9
}

double root(double num)
{
	double root; 
	root=1; 
	while ( fabs(root- num/root ) / root > EPSILON) //1.0E-8 可依精確度需求更改 
	root =(root + num/root)/2; 
	return root;
}


/*--------------------------------------------------------------------------
 * The bisection method.
 *     x, y: the initial value of the root.
 */
double  newton(double *x, double *y)
{
  int     i=0;
  double  err, xnew, ynew, xold, yold;
  double  h, k;
  double  fdx, fdy, gdx, gdy, Delta;
  double  fold, gold;
  double err2;
  FILE *fp;
  //Set up initial conditions.
  xold = *x;
  yold = *y;

  fold = f(xold, yold);
  gold = g(xold, yold);
  fdx = fx(xold, yold);
  fdy = fy(xold, yold);
  gdx = gx(xold, yold);
  gdy = gy(xold, yold);

  Delta = fdx*gdy - fdy*gdx;//fx(xold, yold)*gy(xold, yold)-fy(xold, yold)*gx(xold, yold)
  fp = fopen("data.txt","w");
  xnew = xold ;//修改過 -(
  ynew = yold ;//修改過
  err = 1; // ??error 不知要不要用y的
  fprintf(stderr, "     i            xn			yn               error\n");//修改過 yn
  fprintf(stderr,"------------------------------------------------------------\n");
  //fprintf(fp, "     i            xn			yn               error\n");
  //fprintf(fp,"------------------------------------------------------------\n");
  while(err>=EPSILON){
	  Delta=fabs(fx(xold, yold)*gy(xold, yold)-fy(xold, yold)*gx(xold, yold));
	  if(Delta==0)Delta=EPSILON*10;
   if(i==0)
   {
	   fprintf(stderr,"     %3d \t%lf \t%lf \t---\n", i, xnew, ynew );    //修改過 \t%lf ynew
	   fprintf(fp,"     %3d \t%lf \t%lf \t%lf\n", i, xnew, ynew, err );    //修改過 \t%lf ynew
   }
   else
   {
	   fprintf(stderr,"     %3d \t%lf \t%lf \t%lf\n", i, xnew, ynew , err);
	   fprintf(fp,"     %3d \t%lf \t%lf \t%lf\n", i, xnew, ynew , err);
   }
   xold = xnew; // Save current value.
   yold = ynew; // Save current value.  修改過
   h = -((f(xold, yold)*gy(xold, yold)-g(xold, yold)*fy(xold, yold)))/Delta;
   k = -((-f(xold, yold)*gx(xold, yold)+g(xold, yold)*fx(xold, yold)))/Delta;
   xnew = xold + h; // Compute new values.  修改過xnew = xold - f(xold,yold)/fx(xold,yold);
   ynew = yold + k;
   err = sqrt(fabs((xnew - xold)*(xnew - xold)+(ynew - yold)*(ynew - yold))); // Compute the difference.  (fabs(xnew - xold)+fabs(ynew - yold)); root(fabs((ynew*ynew+xnew*xnew)+(yold*yold+xold*xold)-x*(ynew*yold+xnew*xold))); root(fabs((xnew - xold)*(xnew - xold)+(ynew - yold)*(ynew - yold)))
   xold = xnew; // Save current value.
   yold = ynew; // Save current value.  修改過
   i ++;
  }   
  fprintf(stderr,"------------------------------------------------------------\n"); 
  fprintf(stderr,"     %3d \t%lf \t%lf \t%lf\n", i, xnew, ynew , err);//修改過fprintf(stderr,"     %3d \t%lf \t%lf\n", i, xnew, err);
  //fprintf(fp,"------------------------------------------------------------\n"); 
  fprintf(fp,"     %3d \t%lf \t%lf \t%lf\n", i, xnew, ynew , err);//修改過fprintf(stderr,"     %3d \t%lf \t%lf\n", i, xnew, err);
  fclose(fp);
  return(xnew,ynew);
}


/*----------------------------------------------------------------------------
  * The main procedure.
  */
int main(int argc, char **argv)
{
   double    r;
   double    s;//修改過
   fprintf(stderr," Input the initial value x0=,y0= ");
   fscanf(stdin,"%lf,%lf", &r,&s);
   newton(&r,&s);//修改過 &s增加
   fgetchar();
   system("pause");
}
