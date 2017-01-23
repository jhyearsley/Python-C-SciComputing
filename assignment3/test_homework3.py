import unittest
import numpy

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


#
# Implement your tests here.
#
class TestIntegrate(unittest.TestCase):

    # serial
    def null_test(self):
        self.assertTrue(False)


if __name__ == '__main__':
    #
    # Example timing code: these are on trapz
    #
    print '\n===== Running Timings ====='
    x = linspace(-1, 3, 10000000)
    y = sin(exp(x))

    t1 = time_simps_parallel_chunked(y, x, num_threads=1, chunk_size=1000)
    print 'threads = %d, time = %f'%(1, t1)

    t2 = time_simps_parallel_chunked(y, x, num_threads=2,chunk_size=1000)
    print 'threads = %d, time = %f'%(2, t2)

    t4 = time_simps_parallel_chunked(y, x, num_threads=4,chunk_size=1000)
    print 'threads = %d, time = %f'%(4, t4)

    t8 = time_simps_parallel_chunked(y, x, num_threads=8,chunk_size=1000)
    print 'threads = %d, time = %f'%(8, t8)



    #
    # Run the test suite you created above
    #
    print '\n===== Running Tests ====='
    unittest.main(verbosity=2)


b = linspace(0,1,100);
fb = sin(b);

