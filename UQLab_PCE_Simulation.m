tic; close all; clear all; clc

%{
Polynomial Chaos Expansion Propagation of Satellite Uncertainty using UQLab

This program propagates the uncertainty distribution of a 
satellite in Earth orbit. A Polynomial Chaos Expansion simulation from the
UQLab framework is used based on the given initial Gaussian distribution.
%}

%% Inputs

t0 = 0; %Initial time in seconds
t = 172447.997; %End time in seconds
m0 = [3.9583*10^4 -1.4667*10^4 0.1035*10^4 1.0583 2.8815 0.0842]; %Initial mean vector
p0 = [24287918.0736715,5665435.69165342,894967.886653782,-260.261998968652,1843.15218222622,25.0611461380351;5665435.69165342,257826685.618099,4584696.24234150,-19879.9603291687,247.507838264477,-643.710040075805;894967.886653782,4584696.24234150,150514.823858886,-347.533001229772,64.0785106860803,-7.14358443006258;-260.261998968652,-19879.9603291687,-347.533001229772,1.53503734807762,-0.00580176692199825,0.0499068841013254;1843.15218222622,247.507838264477,64.0785106860803,-0.00580176692199825,0.140297572090407,0.00226834102064119;25.0611461380351,-643.710040075805,-7.14358443006258,0.0499068841013254,0.00226834102064119,0.00244767655339214]*10^-3; %Initial covariance matrix

%% Calculation

% Define model
uqlab
modelopts.mFile = 'UQLab_PCE_Model';
modelopts.isVectorized = false;
modelopts.Parameters = [t0 t];
model = uq_createModel(modelopts);

% Input parameters to model
IOpts.Marginals(1).Type = 'Gaussian';
IOpts.Marginals(1).Parameters = [m0(1) sqrt(p0(1,1))];
IOpts.Marginals(2).Type = 'Gaussian';
IOpts.Marginals(2).Parameters = [m0(2) sqrt(p0(2,2))];
IOpts.Marginals(3).Type = 'Gaussian';
IOpts.Marginals(3).Parameters = [m0(3) sqrt(p0(3,3))];
IOpts.Marginals(4).Type = 'Gaussian';
IOpts.Marginals(4).Parameters = [m0(4) sqrt(p0(4,4))];
IOpts.Marginals(5).Type = 'Gaussian';
IOpts.Marginals(5).Parameters = [m0(5) sqrt(p0(5,5))];
IOpts.Marginals(6).Type = 'Gaussian';
IOpts.Marginals(6).Parameters = [m0(6) sqrt(p0(6,6))];
IOpts.Copula.Type = 'Gaussian';
copula0 = uq_GaussianCopula(corrcov(p0));
IOpts.Copula.Parameters = copula0.Parameters;
input = uq_createInput(IOpts);

% Create metamodel
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Method = 'quadrature';
metamodel = uq_createModel(MetaOpts);

mean = [metamodel.PCE(1).Moments.Mean, metamodel.PCE(2).Moments.Mean, metamodel.PCE(3).Moments.Mean, metamodel.PCE(4).Moments.Mean, metamodel.PCE(5).Moments.Mean, metamodel.PCE(6).Moments.Mean];

fprintf('UQLab Polynomial Chaos Expansion Simulation mean: x = %.3e, y = %.3e, z = %.3e, xdot = %.3e, ydot = %.3e, zdot = %.3e\n',mean)

toc