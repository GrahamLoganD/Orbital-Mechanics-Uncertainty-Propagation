function Y = UQLab_PCE_Model(X, P)
    %{
    Solves the differential equation for the satellite's motion.
        X
            initial state vector [x in km, y in km , z in km, xdot in km/s, ydot in km/s, zdot in km/s]
        P
            parameters [initial time in seconds, final time in seconds]
    %}
    t0 = P(1);
    t = P(2);
    mu = 398600; % Gravitational constant in km^3/s^2
    r = @(x) sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    [~, y] = ode45(@(ti, yi) [yi(4), yi(5), yi(6), -mu * yi(1) / r(yi)^3, -mu * yi(2) / r(yi)^3, -mu * yi(3) / r(yi)^3].', [t0 t], X);
    Y = y(end, 1:end);
end