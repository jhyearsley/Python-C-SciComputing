/*
  integrate.h
  -----------

  Defines routines for numerically integrating (x_i,fval_i) data. That is, given
  a function f, some points along a domain x = [x_0, x_1, ..., x_{N-1}] and
  function values

  fvals = [f(x_0), f(x_1), ..., f(x_{N-1})]

  numerically approximate the integral of f using these function evaluations.
  These are meant to replecate the work done by `scipy.integrate.trapz` and
  `scipy.integrate.simps`, but written in C and using OpenMP.
*/


/*
trapz_serial

Approximates the integral of a function 'f(x)' via a serial implementation of the trapezoidal
rule.

Parameters
----------
fvals : double*
  An array containing the function values of 'f'.

x : double*
  An array containing the respective domain points 'x_i' that the function 'f(x_i)' is evaluated at.

N : int
  The number of points used in the domain.

Returns
-------
double 
  (Return by value.) The function outputs the approximate value of the integral of 'f(x)'

 */
double trapz_serial(double* fvals, double* x, int N);

/*
trapz_parallel

Approximates the integral of a function 'f(x)' via a parallel implementation of the trapezoidal
rule.

Parameters
----------
fvals : double*
  An array containing the function values of 'f(x)' used for the approximation.

x : double*
  An array containing the respective domain points 'x_i' at which 'f(x_i)' are evaluated at.

N : int
  The number of domain points used in the approximation.

num_threads : int
  The number of threads used in the parallelization of the algorithm.

Returns
-------
double
  (Return by value.) The function outputs the approximate value of the integral of 'f(x)'.
 */
double trapz_parallel(double* fvals, double* x, int N, int num_threads);

/*
time_trapz_parallel

A function to calculate the run time of the trapz_parallel algorithm.

Parameters
----------
fvals : double*                                                                                                 
  An array containing the function values of 'f(x)' used for the approximation.                               

x : double*                                                                                                     
  An array containing the respective domain points 'x_i' at which 'f(x_i)' are evaluated at.                  

N : int                                                                                                         
  The number of domain points used in the approximation.                                                      

num_threads : int                                                                                               
  The number of threads used in the parallelization of the algorithm.                                         

Returns                                                           
-------
double
  (Return by value.) The run time of the trapz_parallel algorithm.
 */
double time_trapz_parallel(double* fvals, double* x, int N, int num_threads);

/*
simps_serial

Approximates the integral of the function 'f(x)' on domain points 'x_i'in serial
via Simpsons method.

Parameters
----------
fvals : double*
  An array containing the function values 'f(x_i)' for the domain points 'x_i'.

x : double*
  An array containing the domain points 'x_i'.

Returns
-------
double
  (Return by value.) The approximate integral of 'f(x)' on the specified domain.
 */
double simps_serial(double* fvals, double* x, int N);

/*
simps_parallel

Approximates the integral of the function 'f(x)' on domain points 'x_i' in parallel
via Simpsons method.

Parameters
----------
fvals : double*
  An array containing the function values 'f(x_i)'.

x : double*
  An array containing the domain points 'x_i'.

N : int
  The number of domain points used in the approximation.

num_threads : int
  The number of parallel threads generated.

Returns
-------
double
  (Return by value.) The approximate integral of 'f(x)' on the specified domain.
 */
double simps_parallel(double* fvals, double* x, int N, int num_threads);

/*
time_simps_parallel

A function to calculate the runtime of the simps_parallel algorithm.

Parameters
----------
fvals : double*
  An array containing the function values 'f(x_i)'

x : double*
  An array containing the domain points 'x_i'.

N : int
  The number of domain points used in the approximation.

num_threads : int 
  The number of parallel threads generated.

Returns
-------
double
  (Return by value.) The run time of the simps_parallel algorithm.
 */
double time_simps_parallel(double* fvals, double* x, int N, int num_threads,
                           int repeat);

/*
simps_parallel_chunked

A function to calculate the approximate integral of a function 'f(x)' on domain points 'x_i' in 
parallel which in addition allows the user to specify the amount of work each thread does (chunk size).

Parameters
----------
fvals : double*
  An array containing the function values 'f(x_i)'.

x : double*
  An array containing the domain points 'x_i'.

N : int
  The number of domain points used in the approximation.

num_threads : int
  The number of parallel threads generated.

chunk_size : int
  The number of iterations each thread will perform.

Returns
-------
double
  (Return by value.) The approximate integral of 'f(x)' on the given domain.
 */
double simps_parallel_chunked(double* fvals, double* x, int N, int num_threads,
                              int chunk_size);

/*
time_simps_parallel_chunked

A function to calculate the run time of the simps_parallel_chunked algorithm.

Parameters
----------
fvals : double*                                                                                                
An array containing the function values 'f(x_i)'.                                                          

x : double*                                                                                                   
An array containing the domain points 'x_i'.                                                                

N : int                                                                                                       
The number of domain points used in the approximation.                                                      

num_threads : int                                                                                             
The number of parallel threads generated.                                                                   

chunk_size : int                                                                                              
The number of iterations each thread will perform.                                                          

Returns                                                                                                      
-------
double
  (Return by value.) The run time of the simps_parallel_chunked algorithm.
 */
double time_simps_parallel_chunked(double* fvals, double* x, int N,
                                   int num_threads, int chunk_size, int repeat);
