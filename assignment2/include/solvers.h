#ifndef __homework2_solver_h
#define __homework2_solver_h

#include "linalg.h"

// matrices L, U, and A are all long arrays of size NxN
// b is an array of length N


/*
  solve_lower_triangular

  Solves the matrix equation 'Lx=b' where 'L' is a lower triangular  matrix.
  'L' is an NxN matrix and 'b' is a vector of length N.

  Parameters:
  -----------
  out : double*
    Storage for the solution of the matrix equation.

  L : double*
    A lower triangular matrix

  b : double*
    A vector 'b' of known values, used to determine the solution of 'Lx=b'

  N : int
    The length of the vector 'b'; The # of rows and columns in 'L'

   Returns:
   --------
   out : double*
     (Output by reference.) The vector 'out' is the solution to the matrix equation 'Lx=b'.

 */
void solve_lower_triangular(double* out, double* L, double* b, int N);

/*
  solve_upper_triangular

  Solves the matrix equation 'Ux=b' where 'U' is an NxN upper triangular matrix,
  and 'b' is a vector of length N.

  Parameters:
  -----------
  out : double*
    Storage for the output solution of the matrix equation.

  U : double*
    An upper triangular square matrix.

  b : double*
    A vector 'b' of known values, used to determine the solution of 'Ux=b'.

  N : int
    The length of the vector 'b'; The # of rows and columns of the matrix 'U'.

  Returns:
  --------
  out : double*
    (Return by reference.) The vector 'out' is the solution to the matrix equation 'Ux=b'.
 */
void solve_upper_triangular(double* out, double* U, double* b, int N);

// jacobi returns the number of iterations by value as an `int`. (it also
// returns the solution vector by reference as `out`)
/*
  jacobi

  Iteratively solves the equation Ax=b via the Jacobi method.

  Parameters:
  -----------
  out : double*
    Storage for the output solution of the matrix equation Ax=b.

  A : double*
    A strictly diagonally dominant NxN matrix. It should be noted that
    A is represented as an array of length N*N in this particular algorithm
    *Note* If A is not strictly diagonally dominant the method will not converge.

  b : double*
    An Nx1 vector represented by an array of length N. The vector 'b' of the matrix equation Ax=b.

  N : int
    The # of rows and columns of A; the # of elements in the vector 'b'.

  epsilon: double
    A parameter used to measure convergence criterion. As epsilon becomes smaller the solution
    is a finer approximation.

  Returns:
  --------
  int
    (Return by value.) The function returns by value the number of iterations performed.

  out : double*
    (Return by reference.) The approximate solution to the equation Ax=b.

 */
int jacobi(double* out, double* A, double* b, int N, double epsilon);

// gauss_seidel returns the number of iterations by value as an `int`. (it also
// returns the solution vector by reference as `out`)
/*
  gauss_seidel

  Iteratively solves the matrix equation Ax=b via the Gauss-Seidel method.

  Parameters:
  -----------
  out: double*
    Storage for the output solution

  A : double*
    An array of length N*N to represent the square matrix A of the equation Ax=b.

  b : double*
    An array of length N to represent the vector b of the equation Ax=b.

  N : int
    The # of rows and columns of A; the # of elements of b.

  epsilon : double
    A parameter to tune convergence criterion. As epsilon becomes smaller the approximate
    solution becomes finer.

  Returns:
  --------
  int
    (Return by value.) The function returns by value the number of iterations performed.

  out : double*
    (Return by reference.) The approximate solution to the equation Ax=b.
 */

int gauss_seidel(double* out, double* A, double* b, int N, double epsilon);

#endif
