#include "linalg.h"
#include "solvers.h"
#include  <math.h>
#include <stdlib.h>
// matrices L, U, and A are all long arrays of size NxN
// b is an array of length N




void solve_lower_triangular(double* out, double* L, double* b, int N)
{
  for (int i=0; i<N; ++i)
    out[i] = 0;

  out[0] = b[0]/L[0];
  for (int i=1; i<N; ++i){
    out[i] = b[i];
    for (int j=0; j<i; ++j){
      out[i] -= out[j]*L[i*N+j];
    }
    out[i] = out[i]/L[i*N+i];
  }
}

void solve_upper_triangular(double* out, double* U, double* b, int N)
{
  for (int i=0; i<N; ++i)
    out[i] = 0;

  out[N-1] = b[N-1]/U[N*N-1];
  for (int i=N-2; i>-1; --i){
    out[i] = b[i];
    for (int j=i+1; j<N; ++j){
      out[i] -= out[j]*U[i*N + j];
    } 
    out[i] = out[i]/U[i*N+i];
  }
}


int jacobi(double* out, double* A, double* b, int N, double epsilon)
{
  double* diagA = (double*) malloc(N*N*sizeof(double));
  double* diagAinv = (double*) malloc(N*N*sizeof(double));
  double* restofA = (double*) malloc(N*N*sizeof(double));
  double* rhsvecdiff = (double*) malloc(N*sizeof(double));
  double* restAx = (double*) malloc(N*sizeof(double));
  double* x0 = (double*) malloc(N*sizeof(double));
  double* xnext = (double*) malloc(N*sizeof(double));
  double error;
  double* diff = (double*) malloc(N*sizeof(double));
  int num_iter=0;
  
  for (int i=0; i<N*N; ++i)
    {
      diagA[i] = 0;
    }
  
  for (int i=0; i<N*N; ++i)
    diagAinv[i] = diagA[i];

  for (int i=0; i<N; ++i)
    diagA[i*N+i] = A[i*N+i];

  for (int i=0; i<N; ++i)
    diagAinv[i*N+i] = 1/diagA[i*N+i];
  
  for (int i=0; i<N*N; ++i)
    restofA[i] = A[i] - diagA[i];

  for (int i=0; i<N; ++i)
    x0[i] = 0;


  mat_vec(restAx, restofA, x0, N, N);
  
  vec_sub(rhsvecdiff, b, restAx, N);

  mat_vec(xnext, diagAinv, rhsvecdiff, N, N);

  vec_sub(diff, xnext, x0, N);
  error = vec_norm(diff, N);

  while (error > epsilon)
    {
      num_iter += 1;
      
      for (int i=0; i<N; ++i)
	x0[i] = xnext[i];

      mat_vec(restAx, restofA, x0, N, N);
      
      vec_sub(rhsvecdiff, b, restAx, N);
      
      mat_vec(xnext, diagAinv, rhsvecdiff, N, N);

      vec_sub(diff,xnext, x0, N);
      error= vec_norm(diff, N);
    }


  for (int i=0; i<N; ++i)
    out[i] = xnext[i];
  
  free(diagA);
  free(diagAinv);
  free(restofA);
  free(rhsvecdiff);
  free(restAx);
  free(x0);
  free(xnew);
  free(diff);


  return num_iter;
} 
  
  


int gauss_seidel(double* out, double* A, double* b, int N, double epsilon)
{
  double* upper_A = (double*) malloc(N*N*sizeof(double));
  double* lower_A = (double*) malloc(N*N*sizeof(double));
  double* xold    = (double*) malloc(N*sizeof(double));
  double* lower_A_xold = (double*) malloc(N*sizeof(double));
  double* rhsvecdiff = (double*) malloc(N*sizeof(double));
  double* xnew = (double*) malloc(N*sizeof(double));
  double* diff = (double*) malloc(N*sizeof(double));
  double error;
  int num_iter=0;  

  for (int i=0; i<N; ++i)
    xold[i] = 0;

  for (int i=0; i<N*N; ++i)
    upper_A[i] = 0;

  for (int i=0; i<N*N; ++i)
    lower_A[i] = 0;

  

  for (int i=0; i<N; ++i){
    for (int j=i; j<N; ++j){
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

  while(error > epsilon){
    num_iter += 1;
    for (int i=0; i<N; ++i)
      xold[i] = xnew[i];
    
    mat_vec(lower_A_xold, lower_A, xold, N, N);
    vec_sub(rhsvecdiff, b, lower_A_xold, N);
    solve_upper_triangular(xnew, upper_A, rhsvecdiff, N);
    vec_sub(diff, xnew, xold, N);
    error = vec_norm(diff, N);
  }

  for (int i=0; i<N; ++i)
    out[i] = xnew[i];

  free(upper_A);
  free(lower_A);
  free(xold);
  free(lower_A_xold);
  free(rhsvecdiff);
  free(xnew);
  free(diff);
  return num_iter;
}
