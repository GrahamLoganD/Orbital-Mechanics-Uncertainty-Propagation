function Y = UQLab_PCE_Model(X, P)
    t0 = P(1);
    t = P(2);
    mu=398600; %Gravitational constant in km^3/s^2
    r = @(x) sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    [~,y] = ode45(@(ti,yi) [yi(4) yi(5) yi(6) -mu*yi(1)/r(yi)^3 -mu*yi(2)/r(yi)^3 -mu*yi(3)/r(yi)^3].',[t0 t],X);
    Y = y(end,1:end);
end