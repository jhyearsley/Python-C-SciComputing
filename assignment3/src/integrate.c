#include <omp.h>

double trapz_serial(double* fvals, double* x, int N)
{
  double ans=0;
  for (int i=0; i<(N-1); ++i)
    ans +=(x[i+1] - x[i])*(0.5)*(fvals[i] + fvals[i+1]);

  return ans;
}


double trapz_parallel(double* fvals, double* x, int N, int num_threads)
{
  omp_set_num_threads(num_threads);
  double ans=0;
#pragma omp parallel for reduction(+:ans)
    for (int i=0; i<(N-1); ++i)
    {
      ans += (x[i+1] - x[i])*(0.5)*(fvals[i] + fvals[i+1]);
    }

  return ans;
}


double time_trapz_parallel(double* fvals, double* x, int N, int num_threads)
{
  double end, start = omp_get_wtime();
  trapz_parallel(fvals, x, N, num_threads);
  end = omp_get_wtime();
  return (end - start);
}


double simps_serial(double* fvals, double* x, int N)
{
  double ans=0;
  
  if (N%2 == 0)
    {
      for (int i=0; i<(N-3); i=i+2)
	ans += (x[i+2] - x[i])*(1./6.)*(fvals[i] + 4.*fvals[i+1] + fvals[i+2]);

      for (int i=(N-2); i<(N-1); ++i)
	ans += (x[i+1] - x[i])*(0.5)*(fvals[i] + fvals[i+1]);
    }
  
  else
    {
      for (int i=0; i<(N-2); i=i+2)
	ans += (x[i+2] - x[i])*(1./6.)*(fvals[i] + 4.*fvals[i+1] + fvals[i+2]);
    }
  return ans;
}


double simps_parallel(double* fvals, double* x, int N, int num_threads)
{

  double ans=0;
  
  if (N%2 == 0)
    {
      omp_set_num_threads(num_threads);
      ans = (x[N-1] - x[N-2])*(0.5)*(fvals[N-2] + fvals[N-1]);
#pragma omp parallel for reduction(+:ans)
      for (int i=0; i<(N-3); i=i+2)
	{
	  ans += (x[i+2] - x[i])*(1./6.)*(fvals[i] + 4.*fvals[i+1] + fvals[i+2]);
	}

    }
  
  else
    {
      omp_set_num_threads(num_threads);
      
#pragma omp parallel for reduction(+:ans)
      for (int i=0; i<(N-2); i=i+2)
	  {
	    ans += (x[i+2] - x[i])*(1./6.)*(fvals[i] + 4.*fvals[i+1] + fvals[i+2]);
	  }
    }
  return ans;
}


double time_simps_parallel(double* fvals, double* x, int N, int num_threads,
                           int repeat)
{
  double end, start = omp_get_wtime();
  for (int i=0; i<repeat; ++i)
    simps_parallel(fvals, x, N, num_threads);
  end = omp_get_wtime();
  return (end - start) / (double)repeat;
}


double simps_parallel_chunked(double* fvals, double* x, int N,
                              int num_threads, int chunk_size)
{
  double ans=0;

  if (N%2 == 0)
    {
      omp_set_num_threads(num_threads);
      ans = (x[N-1] - x[N-2])*(0.5)*(fvals[N-2] + fvals[N-1]);

#pragma omp parallel for schedule(static, chunk_size) reduction(+:ans)
      for (int i=0; i<(N-3); i=i+2)
        {
          ans += (x[i+2] - x[i])*(1./6.)*(fvals[i] + 4.*fvals[i+1] + fvals[i+2]);
        }

    }

  else
    {
      omp_set_num_threads(num_threads);

#pragma omp parallel for schedule(static, chunk_size) reduction(+:ans)
      for (int i=0; i<(N-2); i=i+2)
	{
	  ans += (x[i+2] - x[i])*(1./6.)*(fvals[i] + 4.*fvals[i+1] + fvals[i+2]);
	}
    }
  return ans;
}


double time_simps_parallel_chunked(double* fvals, double* x, int N,
                                   int num_threads, int chunk_size,
                                   int repeat)
{
  double end, start = omp_get_wtime();
  for (int i=0; i<repeat; ++i)
    simps_parallel_chunked(fvals, x, N, num_threads, chunk_size);
  end = omp_get_wtime();
  return (end - start) / (double)repeat;
}
