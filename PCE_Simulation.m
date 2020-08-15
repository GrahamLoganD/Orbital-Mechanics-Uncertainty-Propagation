tic; close all; clear all; clc

%{
Polynomial Chaos Expansion Propagation of Satellite Uncertainty

This program propagates the uncertainty distribution of a 
satellite in Earth orbit. A Polynomial Chaos Expansion simulation is used
based on the given initial Gaussian distribution.
%}

%% Inputs

p = 2;  %Maximum order of the polynomial basis
b = 0;  %Extra PCE sampled input uncertainties

t0 = 0; %Initial time in seconds
t = 172447.997; %End time in seconds
m0 = [3.9583*10^4 -1.4667*10^4 0.1035*10^4 1.0583 2.8815 0.0842]; %Initial mean vector
p0 = [24287918.0736715,5665435.69165342,894967.886653782,-260.261998968652,1843.15218222622,25.0611461380351;5665435.69165342,257826685.618099,4584696.24234150,-19879.9603291687,247.507838264477,-643.710040075805;894967.886653782,4584696.24234150,150514.823858886,-347.533001229772,64.0785106860803,-7.14358443006258;-260.261998968652,-19879.9603291687,-347.533001229772,1.53503734807762,-0.00580176692199825,0.0499068841013254;1843.15218222622,247.507838264477,64.0785106860803,-0.00580176692199825,0.140297572090407,0.00226834102064119;25.0611461380351,-643.710040075805,-7.14358443006258,0.0499068841013254,0.00226834102064119,0.00244767655339214]*10^-3; %Initial covariance matrix

%% Calculation

n = length(m0); %Dimensionality
P = factorial(p + n) / (factorial(p) * factorial(n)) - 1;
N2 = factorial(p + n) / (factorial(p) * factorial(n)) + b; %Number of PC sampled input uncertainties

xi = randn([N2 6]); %Sampled input uncertainties
Sx = cholcov(p0);

Y = zeros(N2,n);
for i = 1:N2
    Y(i,:) = phi(t,m0' + Sx*xi(i,:)',t0); %Output uncertainties
end

%Generate matrix of hermite polynomials
H = zeros(N2,1+P);
H(:,1) = ones(N2,1);

for i = 1:N2
    position = 2;
    for j = 1:p
        for k1 = 0:j
            for k2 = 0:j-k1
                for k3 = 0:j-k1-k2
                    for k4 = 0:j-k1-k2-k3
                        for k5 = 0:j-k1-k2-k3-k4
                            k6 = j-k1-k2-k3-k4-k5;
                            H(i,position) = hermiteH(k1,xi(i,1))*hermiteH(k2,xi(i,2))*hermiteH(k3,xi(i,3))*hermiteH(k4,xi(i,4))*hermiteH(k5,xi(i,5))*hermiteH(k6,xi(i,6));
                            position = position + 1;
                        end
                    end
                end
            end
        end
    end
end

C = (H' * H)^-1 * H' * Y; %PC coefficients

mean = C(1,:);
Csub = C(2:end,:)'; %Submatrix of C
cov = Csub * Csub'; %Covariance

fprintf('Polynomial Chaos Expansion Simulation mean: x = %.3e, y = %.3e, z = %.3e, xdot = %.3e, ydot = %.3e, zdot = %.3e\n',mean)

toc

%% Function

function x = phi(t,x0,t0)
    mu = 398600; %Gravitational constant in km^3/s^2
    r = @(x) sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    [~,y] = ode45(@(ti,yi) [yi(4) yi(5) yi(6) -mu*yi(1)/r(yi)^3 -mu*yi(2)/r(yi)^3 -mu*yi(3)/r(yi)^3].',[t0 t],x0);
    x=y(end,1:end);
end