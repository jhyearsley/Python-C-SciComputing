
#ifndef __homework2_linalg_h
#define __homework2_linalg_h

/*
  vec_add

  Computes the sum of two vectors.

  Parameters
  ----------
  out : double*
    Storage for the resulting sum vector.
  v : double*
  w : double*
    The two vectors to sum.
  N : int
    The length of the vectors, `out`, `v`, and `w`.

  Returns
  -------
  out : double*
    (Output by reference.) The sum of `v` and `w`.
*/
void vec_add(double* out, double* v, double* w, int N);


/*
  vec_sub

  Computes the difference of two vectors.

  Parameters
  ----------
  out : double*
    Storage for the resulting difference vector.
  v : double*
  w : double*
    The two vectors used to genereate the difference vector
    (*Note that 'v-w' does not equal 'w-v in general).
  N : int
    The length of 'out', 'v', and 'w'.

  Returns
  -------
  out : double*
    (Output by reference.) The difference of 'v' and 'w'.
 */

void vec_sub(double* out, double* v, double* w, int N);


/*
  vec_norm

  Computes the Euclidean 2-norm of a vector.

  Parameters
  ----------
  v : double*
    The vector you want to compute the norm of.
  N : int
    The length of the vector 'v'.

  Returns
  -------
  vec_norm : double
    (Output by value.) The norm of the vector 'v'.


 */
double vec_norm(double* v, int N);


// represent out, A, and B by arrays of length M*N

/*
  mat_add

  Computes the sum of two matrices.

  Parameters
  ----------
  out : double*
    Storage for the resulting sum matrix.
  A : double*
  B : double*
    The two matrices used to compute the matrix sum.
  N : int
    The length of both matrices. If 'A' and 'B' are M*K matrices
    then 'N=M*K'.

  Returns
  -------
  out : double*
    (Output by reference.) The sum of the matrices 'A' and 'B'.


 */

void mat_add(double* out, double* A, double* B, int M, int N);


// represent A by an array of length M*N

/*
 mat_vec

 Computes the product of a matrix and a vector.

 Parameters
 ----------
 out : double*
   Storage for the resulting matrix vector product.
 A : double*
   The matrix.
 x : double*
   The vector.
 M : int
   The number of rows of the matrix 'A'.
 N : int
   The number of columns of the matrix 'A'; Also the number
   of rows of the vector 'x'.

 Returns
 -------
 out : double*
   (Output by reference.) The product of the matrix 'A' and the vector 'x'.
 

 */
void mat_vec(double* out, double* A, double* x, int M, int N);


// A is MxN, B is NxK, and out is MxK (all as long arrays)

/*
  mat_mat
  
  Computes the product of two matrices

  Parameters
  ----------
  out : double*
    Storage for the resulting matrix product.
  A : double*
  B : double*
    The two matrices to be multiplied.
  M : int
    The number of rows in the matrix 'A'.
  N : int
    The number of columns in the matrix 'A' and the number of rows
    in the matrix 'B'.
  K : int
    The number of columns in the matrix 'B'.

  Returns
  -------
  out : double*
    (Output by reference.) The product of the matrices 'A' and 'B'.

 */
void mat_mat(double* out, double* A, double* B, int M, int N, int K);


#endif
