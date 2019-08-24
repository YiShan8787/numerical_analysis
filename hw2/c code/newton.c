#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

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
  for(i=n-1;i>=0;i--){
    sum = sum*(t-x[i]) + c[i];
	//fprintf(stderr,"i=%d, c[i]= %lf, x[i]= %lf, sum= %lf\n", i, c[i],x[i], sum);
  }
  return sum;
}


/*----------------------------------------------------------
 * The main procedure
 */
void main(int argc, char **argv)
{
  int    n=5;
  double f[]={2.0, 5.0, 6.0, 4.0,2.0};//{0.0, 0.0, 10.0, 10.0, 0.0};
  double x[]={0.0, 5.0, 10.0, 15.0, 0.0};//{0.0, 10.0, 10.0, 0.0, 0.0};
  double c[4];
  double Xw[50],Yw[50];
  double Cx[50],Cy[50],Cw[50];
  double Qt[1000],Qx[1000],Qy[1000],Qw[1000];
  double Px[1000],Py[1000],Pw[50];
  double w[50]={1.0,2.0,1.0,1.0,1.0};
  int     i;
  double  p[4];
  double  T[50]={0.0,sqrt(5.0*5+3.0*3)/(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5)+sqrt(2.0*2+5.0*5)+sqrt(2.0*2+15.0*15)),(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5))/(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5)+sqrt(2.0*2+5.0*5)+sqrt(2.0*2+15.0*15)),(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5)+sqrt(2.0*2+5.0*5))/(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5)+sqrt(2.0*2+5.0*5)+sqrt(2.0*2+15.0*15)),1.0};//squre_chord:{0.0,0.25,0.5,0.75,1}  squre_normal:{0.0,1.0,2.0,3.0,4.0} {0.0,sqrt(5.0*5+3.0*3)/(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5)+sqrt(2.0*2+5.0*5)+sqrt(2.0*2+15.0*15)),(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5))/(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5)+sqrt(2.0*2+5.0*5)+sqrt(2.0*2+15.0*15)),(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5)+sqrt(2.0*2+5.0*5))/(sqrt(5.0*5+3.0*3)+sqrt(1.0+5.0*5)+sqrt(2.0*2+5.0*5)+sqrt(2.0*2+15.0*15)),1.0}
  int num=30;
  FILE *fp;
  srand(time(NULL));
  fp=fopen("data.txt","w");
  //Compute Newton's polynomial.
  for(i=0;i<n;i++)
  {
	  Xw[i]=x[i]*w[i];
	  Yw[i]=f[i]*w[i];
  }
  newton(T, Xw, Cx, n-1);
  newton(T, Yw, Cy, n-1);
  newton(T,w,Cw,n-1);
  //Verify the interpolation conditions
  for(i=0;i<=1000;i++)
  {
	  Qt[i]=(rand()%(int)(T[n-1]-T[0]))+T[0]+(double)rand() / (RAND_MAX+1);//(rand()%(int)(T[n-1]-T[0]))+T[0]+(double)rand() / (RAND_MAX+1);

  }
  for(i=0;i<=1000;i++){
    Qx[i] = Horner(Qt[i], Cx, T, n-1);
	//Py[i] = Horner(Qt[i], Cy, f, n);
    fprintf(stderr,"i=%d, sample= %lf, interpolated= %lf\n", i, Qt[i], Qx[i]);
	//fprintf(stderr,"i=%d, sample= %lf, interpolated= %lf\n", i, Qt[i], Py[i]);

	//fprintf(fp,"%d\t%.10lf\t%.10lf\t%.10lf\n", i, Qt[i], Px[i], Py[i]);
	
  }
  for(i=0;i<=1000;i++){
    //Px[i] = Horner(Qt[i], Cx, x, n);
	Qy[i] = Horner(Qt[i], Cy, T, n-1);
    //fprintf(stderr,"i=%d, sample= %lf, interpolated= %lf\n", i, Qt[i], Px[i]);
	fprintf(stderr,"i=%d, sample= %lf, interpolated= %lf\n", i, Qt[i], Qy[i]);

	//fprintf(fp,"%d\t%.10lf\t%.10lf\t%.10lf\n", i, Qt[i], Px[i], Py[i]);
	
  }
  for(i=0;i<=1000;i++){
    //Px[i] = Horner(Qt[i], Cx, x, n);
	Qw[i] = Horner(Qt[i], Cw, T, n-1);
    //fprintf(stderr,"i=%d, sample= %lf, interpolated= %lf\n", i, Qt[i], Px[i]);
	fprintf(stderr,"i=%d, sample= %lf, interpolated= %lf\n", i, Qt[i], Qw[i]);

	//fprintf(fp,"%d\t%.10lf\t%.10lf\t%.10lf\n", i, Qt[i], Px[i], Py[i]);
	
  }
  for(i=0;i<=1000;i++)
  {
	  Px[i]=Qx[i]/Qw[i];
	  Py[i]=Qy[i]/Qw[i];
	  fprintf(fp,"%d\t%.10lf\t%.10lf\t%.10lf\n", i, Qt[i], Px[i], Py[i]);
  }

  fclose(fp);
}
