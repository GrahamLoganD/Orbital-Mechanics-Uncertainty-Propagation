import math
import time
import numpy
import scipy.stats
import scipy.integrate
import pandas
import seaborn
import matplotlib.pyplot


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


def MC_Sim(N, t0, t, m0, p0, plot):
    '''
    Monte-Carlo Propagation of Satellite Uncertainty

    This program propagates the uncertainty distribution of a 
    satellite in Earth orbit. A Monte-Carlo simulation is used
    based on the given initial Gaussian distribution.
    N
        int, number of random samples
    t0
        float, initial time in seconds
    t
        float, end time in seconds
    m0
        list of floats, initial mean of state vector distribution [x0 in km, y0 in km, z0 in km, xdot0 in km/s, ydot0 in km/s, zdot0 in km/s]
    p0
        numpy matrix, initial covariance matrix of state vector distribution
    plot
        bool, if true a plot of the final distribution will be generated
    '''

    n = len(p0) # Dimensionality

    x0 = scipy.stats.multivariate_normal.rvs(mean=m0, cov=p0, size=N) # Random samples from initial distribution

    # Solve for final states
    x = numpy.zeros([N, n])
    for i in numpy.arange(0, N):
        x[i, :] = phi(t, x0[i, :], t0)

    m = numpy.mean(x, axis=0) # Mean of final distribution

    # Generate plot of final distribution
    if plot:
        df = pandas.DataFrame(x, columns=[r"$x$ (km)", r"$y$ (km)", r"$z$ (km)", r"$\dotx$ (km/s)", r"$\doty$ (km/s)", r"$\dotz$ (km/s)"])
        seaborn.set()
        seaborn.jointplot(x=r"$x$ (km)", y=r"$y$ (km)", data=df, kind="kde")
        seaborn.jointplot(x=r"$x$ (km)", y=r"$z$ (km)", data=df, kind="kde")
        seaborn.jointplot(x=r"$\dotx$ (km/s)", y=r"$\doty$ (km/s)", data=df, kind="kde")
        seaborn.jointplot(x=r"$\dotx$ (km/s)", y=r"$\dotz$ (km/s)", data=df, kind="kde")
        matplotlib.pyplot.show()

    return m


if __name__ == "__main__":
    start = time.time()
    N = 50000
    t0 = 0
    t = 172447.997
    m0 = numpy.array([3.9583 * (10 ** 4), -1.4667 * (10 ** 4), 0.1035 * (10 ** 4), 1.0583, 2.8815, 0.0842])
    p0 = numpy.array([[ 24287918.0736715,  5665435.69165342,  894967.886653782,    -260.261998968652,     1843.15218222622,    25.0611461380351],
                    [ 5665435.69165342,  257826685.618099,  4584696.24234150,    -19879.9603291687,     247.507838264477,   -643.710040075805],
                    [ 894967.886653782,  4584696.24234150,  150514.823858886,    -347.533001229772,     64.0785106860803,   -7.14358443006258],
                    [-260.261998968652, -19879.9603291687, -347.533001229772,     1.53503734807762, -0.00580176692199825,  0.0499068841013254],
                    [ 1843.15218222622,  247.507838264477,  64.0785106860803, -0.00580176692199825,    0.140297572090407, 0.00226834102064119],
                    [ 25.0611461380351, -643.710040075805, -7.14358443006258,   0.0499068841013254,  0.00226834102064119, 0.00244767655339214]]) * (10 ** -3)
    m = MC_Sim(N, t0, t, m0, p0, True)          
    end = time.time()
    print('Monte Carlo Simulation mean: x = {:.4e}, y = {:.4e}, z = {:.4e}, xdot = {:.4e}, ydot = {:.4e}, zdot = {:.4e}'.format(*m.tolist()))
    print('Elapsed time is {} seconds.'.format(end-start))