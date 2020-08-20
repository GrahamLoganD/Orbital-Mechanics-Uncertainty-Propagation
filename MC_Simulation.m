tic; close all; clear all; clc

%{
Monte-Carlo Propagation of Satellite Uncertainty

This program propagates the uncertainty distribution of a 
satellite in Earth orbit. A Monte-Carlo simulation is used
based on the given initial Gaussian distribution.
%}

%% Inputs
N = 50000; % Number of random samples

t0 = 0; % Initial time in seconds
t = 172447.997; % End time in seconds
m0 = [3.9583*10^4 -1.4667*10^4 0.1035*10^4 1.0583 2.8815 0.0842]; % Initial mean of state vector distribution [x0 in km, y0 in km, z0 in km, xdot0 in km/s, ydot0 in km/s, zdot0 in km/s]
p0 = [ 24287918.0736715,  5665435.69165342,  894967.886653782, -260.261998968652,  1843.15218222622,  25.0611461380351;
       5665435.69165342,  257826685.618099,  4584696.24234150, -19879.9603291687,  247.507838264477, -643.710040075805;
       894967.886653782,  4584696.24234150,  150514.823858886, -347.533001229772,  64.0785106860803, -7.14358443006258;
      -260.261998968652, -19879.9603291687, -347.533001229772,  1.53503734807762, -0.00580176692199,  0.04990688410132;
       1843.15218222622,  247.507838264477,  64.0785106860803, -0.00580176692199,  0.14029757209040,  0.00226834102064;
       25.0611461380351, -643.710040075805, -7.14358443006258,  0.04990688410132,  0.00226834102064,  0.00244767655339]*10^-3; % Initial covariance matrix of state vector distribution

%% Dimensionality
n = length(m0);

%% Generate random samples from initial distribution
x0 = mvnrnd(m0, p0, N);

%% Solve for final states
x = zeros(N, n);
for i = 1:N
    x(i, :) = phi(t, x0(i, :), t0);
end

%% Output
fprintf('Monte Carlo Simulation mean: x = %.3e, y = %.3e, z = %.3e, xdot = %.3e, ydot = %.3e, zdot = %.3e\n', mean(x, 1))

toc

%% Plotting

% Plot y versus x
figure('Name', 'f1')
s1 = scatterhist(x(:, 1), x(:, 2), 'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1, 'Marker', '.none');
xlabel('x (km)')
ylabel('y (km)')

% Plot z versus x
figure('Name', 'f2')
s2 = scatterhist(x(:, 1), x(:, 3), 'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1, 'Marker', '.none');
xlabel('x (km)')
ylabel('z (km)')

% Plot ydot versus xdot
figure('Name', 'f3')
s3 = scatterhist(x(:, 4), x(:, 5), 'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1, 'Marker', '.none');
xlabel('$\dot{x}$ (km/s)', 'Interpreter', 'latex')
ylabel('$\dot{y}$ (km/s)', 'Interpreter', 'latex')

% Plot zdot versus xdot
figure('Name', 'f4')
s4 = scatterhist(x(:, 4), x(:, 6), 'Location', 'SouthEast', 'Direction', 'out', 'Kernel', 'on', 'LineWidth', 1, 'MarkerSize', 1, 'Marker', '.none');
xlabel('$\dot{x}$ (km/s)', 'Interpreter', 'latex')
ylabel('$\dot{z}$ (km/s)', 'Interpreter', 'latex')

% Combine plots
figure('Name', 'f5')
u1 = uipanel('position', [0 .5 .5 .5]);
set(s1, 'parent', u1)
u2 = uipanel('position', [.5 .5 .5 .5]);
set(s2, 'parent', u2)
u3 = uipanel('position', [0 0 .5 .5]);
set(s3, 'parent', u3)
u4 = uipanel('position', [.5 0 .5 .5]);
set(s4, 'parent', u4)
close f1 f2 f3 f4

%% Modeling function

function x = phi(t, x0, t0)
    %{
    Solves the differential equation for the satellite's motion.
        t
            final time in seconds
        x0
            initial state vector [x0 in km, y0 in km , z0 in km, xdot0 in km/s, ydot0 in km/s, zdot0 in km/s]
        t0
            initial time in seconds
    %}

    mu = 398600; % Gravitational constant in km^3/s^2
    
    r = @(x) sqrt(x(1)^2 + x(2)^2 + x(3)^2); % Function for satellite's distance
    
    [~, y] = ode45(@(ti, yi) [yi(4), yi(5), yi(6), -mu * yi(1) / r(yi)^3, -mu * yi(2) / r(yi)^3, -mu * yi(3) / r(yi)^3].', [t0 t], x0);
    x = y(end, 1:end);
end