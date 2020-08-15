tic; close all; clear all; clc

%{
Monte-Carlo Propagation of Satellite Uncertainty

This program propagates the uncertainty distribution of a 
satellite in Earth orbit. A Monte-Carlo simulation is used
based on the given initial Gaussian distribution.
%}

%% Inputs

N = 50000; %Number of random samples

t0 = 0; %Initial time in seconds
t = 172447.997; %End time in seconds
m0 = [3.9583*10^4 -1.4667*10^4 0.1035*10^4 1.0583 2.8815 0.0842]; %Initial mean vector
p0 = [24287918.0736715,5665435.69165342,894967.886653782,-260.261998968652,1843.15218222622,25.0611461380351;5665435.69165342,257826685.618099,4584696.24234150,-19879.9603291687,247.507838264477,-643.710040075805;894967.886653782,4584696.24234150,150514.823858886,-347.533001229772,64.0785106860803,-7.14358443006258;-260.261998968652,-19879.9603291687,-347.533001229772,1.53503734807762,-0.00580176692199825,0.0499068841013254;1843.15218222622,247.507838264477,64.0785106860803,-0.00580176692199825,0.140297572090407,0.00226834102064119;25.0611461380351,-643.710040075805,-7.14358443006258,0.0499068841013254,0.00226834102064119,0.00244767655339214]*10^-3; %Initial covariance matrix

%% Calculation

n = length(m0); %Dimensionality

x0 = mvnrnd(m0,p0,N);
x = zeros(N,n);

for i = 1:N
    x(i,:) = phi(t,x0(i,:),t0);
end

fprintf('Monte Carlo Simulation mean: x = %.3e, y = %.3e, z = %.3e, xdot = %.3e, ydot = %.3e, zdot = %.3e\n',mean(x,1))

toc

%% Plotting

figure('Name','f1')
s1=scatterhist(x(:,1),x(:,2),'Location','SouthEast','Direction','out','Kernel','on','LineWidth',1,'MarkerSize',1,'Marker','.none');
xlabel('x (km)')
ylabel('y (km)')

figure('Name','f2')
s2=scatterhist(x(:,1),x(:,3),'Location','SouthEast','Direction','out','Kernel','on','LineWidth',1,'MarkerSize',1,'Marker','.none');
xlabel('x (km)')
ylabel('z (km)')

figure('Name','f3')
s3=scatterhist(x(:,4),x(:,5),'Location','SouthEast','Direction','out','Kernel','on','LineWidth',1,'MarkerSize',1,'Marker','.none');
xlabel('$\dot{x}$ (km/s)','Interpreter','latex')
ylabel('$\dot{y}$ (km/s)','Interpreter','latex')

figure('Name','f4')
s4=scatterhist(x(:,4),x(:,6),'Location','SouthEast','Direction','out','Kernel','on','LineWidth',1,'MarkerSize',1,'Marker','.none');
xlabel('$\dot{x}$ (km/s)','Interpreter','latex')
ylabel('$\dot{z}$ (km/s)','Interpreter','latex')

figure('Name','f5')
u1 = uipanel('position',[0 .5 .5 .5]);
set(s1,'parent',u1)
u2 = uipanel('position',[.5 .5 .5 .5]);
set(s2,'parent',u2)
u3 = uipanel('position',[0 0 .5 .5]);
set(s3,'parent',u3)
u4 = uipanel('position',[.5 0 .5 .5]);
set(s4,'parent',u4)
close f1 f2 f3 f4

%% Function

function x = phi(t,x0,t0)
    mu = 398600; %Gravitational constant in km^3/s^2
    r = @(x) sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    [~,y] = ode45(@(ti,yi) [yi(4) yi(5) yi(6) -mu*yi(1)/r(yi)^3 -mu*yi(2)/r(yi)^3 -mu*yi(3)/r(yi)^3].',[t0 t],x0);
    x=y(end,1:end);
end