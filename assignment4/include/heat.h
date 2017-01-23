#ifndef __homework_4_heat_h
#define __homework_4_heat_h

#include <mpi.h>

/*
  heat_serial

  A serial version of solving the periodic heat equation on a uniform grid. Use
  for comparison with the MPI distributed version.

  The forward Euler method for numerical solving PDEs is unstable for large
  (even only moderately) values of `dt`. Make sure `dt` is small. Examples are
  given in the Python code that calls this function.

  Parameters
  ----------
  u : double*
      Input data array. Represents the heat / temperature across a 1d rod.
      Assumes that the data points are equally spaced in the space domain.
  Nx : size_t
      The length of `u`.
  dt : double
      The size of the time step to take with each iteration.
  Nt : size_t
      Number of time-steps to perform.

  Returns
  -------
  None
      The array `u` is modified in place.

*/
void heat_serial(double* u, double dx, size_t Nx, double dt, size_t Nt);

/*
  heat_parallel

  A parallel version of solving the periodic heat equation on uniform grid.
  The algorithm makes use of openMPI to distribute the domain to be solved
  by a specified amount of processes which each make use of the forward
  Euler method.

  *NOTE* The forward Euler method for numerically solving PDEs is unstable 
  for large (even only moderately) values of `dt`. Make sure `dt` is small.

  Parameters
  ----------
  uk : double*
    Input data representing the local heat/temp across a 1-d rod; the local
    chunk of values used by each process in the integration scheme.

  dx : double
    The discretization size of the domain.

  Nx : size_t
    The length of `uk`.

  dt : double
    The size of the time step to take with each iteration.

  Nt : size_t
    The number of iterations to perform in the time-stepping routine.

  comm : MPI_Comm
    The communicator number that corresponds to the process working on the
    chunk `uk`.

  Returns
  -------
  none
    The array `uk` is modified in place.
 */
void heat_parallel(double* uk, double dx, size_t Nx, double dt, size_t Nt,
                   MPI_Comm comm);

#endif
