import math
import time
import numpy
import scipy.integrate
import scipy.linalg
import scipy.stats


def phi(t, x0, t0):
    '''
        Solves the differential equation for the satellite's motion.
        t
            float, final time in seconds
        x0
            list of floats, initial state vector [x in km, y in km , z in km, xdot in km/s, ydot in km/s, zdot in km/s]
        t0
            float, initial time in seconds
    '''

    mu = 398600 # Gravitational constant in km^3/s^2

    r = lambda x: math.sqrt((x[0] ** 2) + (x[1] ** 2) + (x[2] ** 2)) # Function for satellite's distance

    sol = scipy.integrate.odeint(func=lambda t, y: numpy.array([y[3], y[4], y[5], (-mu * y[0]) / (r(y) ** 3), (-mu * y[1]) / (r(y) ** 3), (-mu * y[2]) / (r(y) ** 3)]), t=[t0, t], y0=x0, tfirst=True)

    return sol[-1, :]


def PCE_Sim(p, b, t0, t, m0, p0):
    '''
    Polynomial Chaos Expansion Propagation of Satellite Uncertainty

    This program propagates the uncertainty distribution of a satellite in 
    Earth orbit. A Polynomial Chaos Expansion simulation is used based on the 
    given initial Gaussian distribution. The method is non-intrusive with 
    least-squares regression.
    p
        int, maximum order of the polynomial basis
    b
        int, extra PCE sampled input uncertainties
    t0
        float, initial time in seconds
    t
        float, end time in seconds
    m0
        list of floats, initial mean of state vector distribution [x0 in km, y0 in km, z0 in km, xdot0 in km/s, ydot0 in km/s, zdot0 in km/s]
    p0
        numpy matrix, initial covariance matrix of state vector distribution
    '''

    # Calculate constants
    n = len(p0) # Dimensionality
    P = int(math.factorial(p + n) / (math.factorial(p) * math.factorial(n)) - 1) # P + 1 is the width of the hermitian matrix
    N2 = int(math.factorial(p + n) / (math.factorial(p) * math.factorial(n)) + b) # Number of PC sampled input uncertainties

    xi = numpy.array(scipy.stats.multivariate_normal.rvs(size=(N2, n))) # Sampled input uncertainties

    # Solve for output uncertanties
    Sx = numpy.matrix(scipy.linalg.cholesky(p0)) # Cholesky decomposition of the initial covariance matrix
    Y = numpy.zeros((N2, n))
    for i in numpy.arange(1, N2):
        Y[i, :] = phi(t, numpy.ravel(numpy.transpose([m0]) + (Sx * numpy.transpose([xi[i, :]]))), t0) # Output uncertainties

    # Generate matrix of hermite polynomials
    H = numpy.matrix(numpy.zeros((N2, 1 + P)))
    H[:, 0] = 1
    for i in numpy.arange(0, N2):
        position = 1
        for j in numpy.arange(1, p + 1):
            for k1 in numpy.arange(0, j + 1):
                for k2 in numpy.arange(0, j - k1 + 1):
                    for k3 in numpy.arange(0, j - k1 - k2 + 1):
                        for k4 in numpy.arange(0, j - k1 - k2 - k3 + 1):
                            for k5 in numpy.arange(0, j - k1 - k2 - k3 - k4 + 1):
                                k6 = j - k1 - k2 - k3 - k4 - k5
                                H[i, position] = scipy.special.eval_hermite(k1, xi[i, 0]) * scipy.special.eval_hermite(k2, xi[i, 1]) * scipy.special.eval_hermite(k3, xi[i, 2]) * scipy.special.eval_hermite(k4, xi[i, 3]) * scipy.special.eval_hermite(k5, xi[i, 4]) * scipy.special.eval_hermite(k6, xi[i, 5])
                                position = position + 1

    C = numpy.linalg.inv(H.T * H) * H.T * Y # Least-squares calculation of PC coefficients

    mean = numpy.ravel(C[0, :]) # Mean of final distribution
    return mean


if __name__ == "__main__":
    start = time.time()
    p = 2
    b = 0
    t0 = 0
    t = 172447.997
    m0 = [3.9583 * (10 ** 4), -1.4667 * (10 ** 4), 0.1035 * (10 ** 4), 1.0583, 2.8815, 0.0842]
    p0 = numpy.matrix([[ 24287918.0736715,  5665435.69165342,  894967.886653782,    -260.261998968652,     1843.15218222622,    25.0611461380351],
                       [ 5665435.69165342,  257826685.618099,  4584696.24234150,    -19879.9603291687,     247.507838264477,   -643.710040075805],
                       [ 894967.886653782,  4584696.24234150,  150514.823858886,    -347.533001229772,     64.0785106860803,   -7.14358443006258],
                       [-260.261998968652, -19879.9603291687, -347.533001229772,     1.53503734807762, -0.00580176692199825,  0.0499068841013254],
                       [ 1843.15218222622,  247.507838264477,  64.0785106860803, -0.00580176692199825,    0.140297572090407, 0.00226834102064119],
                       [ 25.0611461380351, -643.710040075805, -7.14358443006258,   0.0499068841013254,  0.00226834102064119, 0.00244767655339214]]) * (10 ** -3)
    mean = PCE_Sim(p, b, t0, t, m0, p0)
    end = time.time()
    print('Polynomial Chaos Expansion Simulation mean: x = {:.4e}, y = {:.4e}, z = {:.4e}, xdot = {:.4e}, ydot = {:.4e}, zdot = {:.4e}'.format(*mean.tolist()))
    print('Elapsed time is {} seconds.'.format(end-start))