/*
  example.c

  An example file for directly interacting with your C library. Run `make lib`
  to create the dynamic library. From the top level directory (the directory
  above this one named "homework2-githubusername") run

  $ gcc ctests/example.c -Llib -lhomework2 -Iinclude -Wl,-rpath,./lib -lm -o example
  $ ./example

  (See Issue #7 in homework2. Thank you to @rachka, @jlombs, and @shinwookang
  for debugging.)

*/

#include <stdio.h>
#include "linalg.h"
#include "solvers.h"




int main(int argc, char** argv)
{
  double U[4] = {8., 2., 0., 3.};
  double B[4];
  double b[2] = {3., 2.};
  double out1[2];
  double out2[2];
  int num_iter=0;
  double epsilon = 0.00000000001;;
  int N=2;
  double diagA[4];
  double upper_A[4];
  double lower_A[4];
  double xold[2];
  double lower_A_xold[2];
  double rhsvecdiff[2];
  double diff[2];
  double xnew[2];
  double error;
  double diagAinv[4];
  double out[2];
  /*
  num_iter = gauss_seidel(out1, A, b, N, epsilon);
  
  for (int i=0; i<N; ++i)
    xnew[i] = 0;
  

  for (int i=0; i<N*N; ++i)
    upper_A[i] = 0;

  printf("\nGauss ish\n");
  for (int i=0; i<N; ++i)
    printf("%f\n",out1[i]);

  
  printf("N = %d\n",N);
  for (int i=0; i<N; ++i)
    xold[i] = 0;

  for (int i=0; i<N; ++i){
    printf("i = %d\n",i);
    for (int j=i; j<N; ++j){
      printf("j = %d\n",j);
      upper_A[i*N + j] = A[i*N + j];
    }
  }

  for (int i=0; i<N*N; ++i)
    lower_A[i] = A[i] - upper_A[i];

  mat_vec(lower_A_xold, lower_A, xold, N, N);
  vec_sub(rhsvecdiff, b, lower_A_xold, N);
  solve_upper_triangular(xnew, upper_A, rhsvecdiff, N);
  vec_sub(diff, xnew, xold, N);
  error = vec_norm(diff, N);
  
  printf("\n\n Upper A = \n");
  for (int i=0; i<N*N; ++i)
    printf("%f\n", upper_A[i]);


  printf("\nLower_A = \n");
  for (int i=0; i<N*N; ++i)
    printf("%f\n",lower_A[i]);
  printf("\n%d\n",num_iter);
    

  printf("\n X_old = \n");
  for (int i=0; i<N; ++i)
    printf("%f\n",xold[i]);

  printf("\n Lower_A_xold = \n");
  for (int i=0; i<N; ++i)
    printf("%f\n", lower_A_xold[i]);

  printf("\nRhsvecdiff = \n");
  for (int i=0; i<N; ++i)
    printf("%f\n",rhsvecdiff[i]);

  printf("\n Xnew = \n");
  for (int i=0; i<N; ++i)
    printf("%f\n",xnew[i]);
  */

  printf("\n\nNeed space for debugging upper solver!\n");
  for (int i=0; i<N; ++i)
    out[i] = 0;

  out[N-1] = b[N-1]/U[N*N-1];
  printf("out[N-1] = %f\n",out[N-1]);
  for (int i=N-2; i>-1; --i){
    out[i] = b[i];
    for (int j=i+1; j<N; ++j){
      out[i] -= out[j]*U[i*N + j];
    }
    out[i] = out[i]/U[i*N+i];
  }
}
 
