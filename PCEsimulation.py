'''
Polynomial Chaos Expansion Propagation of Satellite Uncertainty

This program propagates the uncertainty distribution of a 
satellite in Earth orbit. A polynomial chaos expansion simulation is used
based on the given initial Gaussian distribution.
'''

import math
import time
import numpy
import scipy.stats
import scipy.linalg
start = time.time()

def phi(t,x0,t0):
    import math
    import numpy
    from scipy.integrate import odeint

    mu = 398600 #Gravitational constant in km^3/s^2
    r = lambda x: math.sqrt((x[0] ** 2) + (x[1] ** 2) + (x[2] ** 2))

    sol = odeint(func = lambda t,y: numpy.array([y[3], y[4], y[5], (-mu * y[0]) / (r(y) ** 3), (-mu * y[1]) / (r(y) ** 3), (-mu * y[2]) / (r(y) ** 3)]), t = [t0, t], y0 = x0, tfirst = True)

    return sol[-1,:]

#   INPUTS:

p = 2  #Maximum order of the polynomial basis
b = 0  #Extra PCE sampled input uncertainties

t0 = 0 #Initial time in seconds
t = 172447.997 #End time in seconds
m0 = numpy.array([3.9583 * (10 ** 4), -1.4667 * (10 ** 4), 0.1035 * (10 ** 4), 1.0583, 2.8815, 0.0842]) #Initial mean vector
p0 = numpy.matrix([[ 24287918.0736715,  5665435.69165342,  894967.886653782,    -260.261998968652,     1843.15218222622,    25.0611461380351],
                  [ 5665435.69165342,  257826685.618099,  4584696.24234150,    -19879.9603291687,     247.507838264477,   -643.710040075805],
                  [ 894967.886653782,  4584696.24234150,  150514.823858886,    -347.533001229772,     64.0785106860803,   -7.14358443006258],
                  [-260.261998968652, -19879.9603291687, -347.533001229772,     1.53503734807762, -0.00580176692199825,  0.0499068841013254],
                  [ 1843.15218222622,  247.507838264477,  64.0785106860803, -0.00580176692199825,    0.140297572090407, 0.00226834102064119],
                  [ 25.0611461380351, -643.710040075805, -7.14358443006258,   0.0499068841013254,  0.00226834102064119, 0.00244767655339214]]) * (10 ** -3) #Initial covariance matrix


#   CALCULATION:

n = len(p0) #Dimensionality
P = int(math.factorial(p + n) / (math.factorial(p) * math.factorial(n)) - 1)
N2 = int(math.factorial(p + n) / (math.factorial(p) * math.factorial(n)) + b) #Number of PC sampled input uncertainties

xi = numpy.array(scipy.stats.multivariate_normal.rvs(size = (N2,6))) #Sampled input uncertainties
Sx = numpy.matrix(scipy.linalg.cholesky(p0))

Y = numpy.zeros((N2,n))
for i in numpy.arange(1,N2):
    Y[i,:] = phi(t,numpy.ravel(numpy.transpose([m0]) + (Sx * numpy.transpose([xi[i,:]]))),t0) #Output uncertainties

#Generate matrix of hermite polynomials
H = numpy.matrix(numpy.zeros((N2,1+P)))
H[:,0] = 1

for i in numpy.arange(0,N2):
    position = 1
    for j in numpy.arange(1,p + 1):
        for k1 in numpy.arange(0,j + 1):
            for k2 in numpy.arange(0,j - k1 + 1):
                for k3 in numpy.arange(0,j - k1 - k2 + 1):
                    for k4 in numpy.arange(0,j - k1 - k2 - k3 + 1):
                        for k5 in numpy.arange(0,j - k1 - k2 - k3 - k4 + 1):
                            k6 = j - k1 - k2 - k3 - k4 - k5
                            H[i,position] = scipy.special.eval_hermite(k1,xi[i,0]) * scipy.special.eval_hermite(k2,xi[i,1]) * scipy.special.eval_hermite(k3,xi[i,2]) * scipy.special.eval_hermite(k4,xi[i,3]) * scipy.special.eval_hermite(k5,xi[i,4]) * scipy.special.eval_hermite(k6,xi[i,5])
                            position = position + 1

C = numpy.linalg.inv(H.T * H) * H.T * Y #PC coefficients

mean = numpy.ravel(C[0,:])
Csub = C[1:-1,:].T
cov = Csub * Csub.T #Covariance

print('Polynomial Chaos Expansion Simulation mean: x = {:.4e}, y = {:.4e}, z = {:.4e}, xdot = {:.4e}, ydot = {:.4e}, zdot = {:.4e}'.format(*mean.tolist()))

end = time.time()
print('Elapsed time is {} seconds.'.format(end-start))