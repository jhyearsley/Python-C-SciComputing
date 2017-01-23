import unittest
import numpy 
import matplotlib.pyplot as plt


from numpy import array, linspace, pi, sin, cos, exp
from scipy.integrate import trapz, simps

from homework3 import (
    trapz_serial,
    trapz_parallel,
    time_trapz_parallel,
    simps_serial,
    simps_parallel,
    time_simps_parallel,
    simps_parallel_chunked,
    time_simps_parallel_chunked,
)


N = 2**20
x = numpy.linspace(-1,3,N)
y = sin(exp(x))

chunk_sizes = [];
timing_nthreads1 = [];
timing_nthreads2 = [];
timing_nthreads4 = [];
timing_nthreads8 = [];
timing_nthreads300 = [];

for i in range(0,20,2):
    chunk_sizes.append(2**i);

for i in chunk_sizes:
    timing_nthreads1.append(time_simps_parallel_chunked(y, x, num_threads=1, repeat=100, chunk_size=i))
    timing_nthreads2.append(time_simps_parallel_chunked(y, x, num_threads=2, repeat=100, chunk_size=i))
    timing_nthreads4.append(time_simps_parallel_chunked(y, x, num_threads=4, repeat=100, chunk_size=i))
    timing_nthreads8.append(time_simps_parallel_chunked(y, x, num_threads=8, repeat=100, chunk_size=i))
    timing_nthreads300.append(time_simps_parallel_chunked(y, x, num_threads=300, repeat=100, chunk_size=i))

plt.semilogx(chunk_sizes, timing_nthreads1, color="k", marker='o', label="1 Thread")
plt.semilogx(chunk_sizes, timing_nthreads2, color="g", marker='o', label="2 Threads")
plt.semilogx(chunk_sizes, timing_nthreads4, color="b", marker='o', label="4 Threads")
plt.semilogx(chunk_sizes, timing_nthreads8, color="r", marker='o', label="8 Threads")
plt.semilogx(chunk_sizes, timing_nthreads300, color="m", marker='o', label="300 Threads")

plt.title("Parallel Simpson's Rule Timings")
plt.xlabel("Chunk Size")
plt.ylabel("Time (sec)")
plt.legend()
plt.show()
