import unittest
import numpy
from numpy import array, ones, diag, arange
from numpy.linalg import norm
from numpy.random import randn
from scipy.linalg import solve_triangular as scipy_solve_triangular
from scipy.linalg import solve as scipy_solve

# import the Python wrappers defined in homework2.wrappers
from homework2 import (
    vec_add,
    vec_sub,
    vec_norm,
    mat_add,
    mat_vec,
    mat_mat,
    solve_lower_triangular,
    solve_upper_triangular,
    jacobi,
    gauss_seidel,
)

# these might be useful for your tests
import numpy
from numpy import array, ones, diag, arange
from numpy.linalg import norm
from numpy.random import randn
from scipy.linalg import solve_triangular as scipy_solve_triangular
from scipy.linalg import solve as scipy_solve

# this also might be useful for your tests
def five_diagonal_system(n):
    """Returns the 5-diagonal matrix A and vector b used in some of the Jacobi and
    Gauss-Seidel tests.

    Parameters
    ----------
    n : int
        The size of the system.

    Returns
    -------
    A : array (matrix)
        The system matrix.
    b : array (vector)
        The system vector.
    """
    v = 5*ones(n)
    w1 = -ones(n-1)
    w2 = -ones(n-2)
    A = diag(v) + diag(w2,k=-2) + diag(w1,k=-1) + diag(w1,k=1) + diag(w2,k=2)
    b = arange(n)
    return A,b


class TestLinalg(unittest.TestCase):
    def test_vec_add(self):
        x = array([1,2,3], dtype=numpy.double)
        y = array([4,5,6], dtype=numpy.double)
        z = vec_add(x,y)
        error = norm(z - (x+y))
        self.assertAlmostEqual(error, 0)

    def test_mat_mat(self):
        A = randn(8,8)  # create two random matrices
        B = randn(8,8)
        C = mat_mat(A,B)  # product using C code (wrapped by homework2.mat_mat)
        C_actual = numpy.dot(A,B)  # product using Numpy
        error = norm(C - C_actual)
        self.assertLess(error, 1e-8)  # account for floating point error

class TestSolver(unittest.TestCase):
    pass # delete this once you write some tests here



def time_gauss_seidel(n, number=1):
    """Returns the amount of time (in seconds) for gauss_seidel to execute for the
    given system size `n`.

    Feel free to use to timing purposes.

    Parameters
    ----------
    n : int
        The size of the test system to solve.
    number : int
        (Optional) The number of times to run the test. Useful for computing an
        average runtime.

    Returns
    -------
    avg_time : float
        The time, in seconds, it took to run gauss_seidel().
    """
    from timeit import timeit

    # write the setup (non-timed but necessary) code as a string
    s = '''
from numpy.random import randn
from homework2 import gauss_seidel
from test_homework2 import five_diagonal_system
N = %d
A,b = five_diagonal_system(N)
'''%(n)
    total_time = timeit('gauss_seidel(A,b)', setup=s, number=number)
    avg_time = total_time / number
    return avg_time

def time_jacobi(n, number=1):
    """Returns the amount of time (in seconds) for gauss_seidel to execute for the
    given system size `n`.

    Feel free to use to timing purposes.

    Parameters
    ----------
    n : int
        The size of the test system to solve.
    number : int
        (Optional) The number of times to run the test. Useful for computing an
        average runtime.

    Returns
    -------
    avg_time : float
        The time, in seconds, it took to run gauss_seidel().
    """
    from timeit import timeit

    # write the setup (non-timed but necessary) code as a string
    s = '''
from numpy.random import randn
from homework2 import jacobi
from test_homework2 import five_diagonal_system
N = %d
A,b = five_diagonal_system(N)
'''%(n)
    total_time = timeit('jacobi(A,b)', setup=s, number=number)
    avg_time = total_time / number
    return avg_time

A,b = five_diagonal_system(800);
j, num = jacobi(A,b);
g, num = gauss_seidel(A,b);
s = scipy_solve(A,b);

print norm(s)
print norm(g)
print norm(j)
if __name__ == '__main__':
    #print '\n===== Timings ====='
    #n = 2**8
    #t = time_gauss_seidel(n)
    #print 'gauss_seidel(%d): %f'%(n, t)

    print '\n===== Running Tests ====='
    unittest.main(verbosity=2)


