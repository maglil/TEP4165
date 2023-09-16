function T = upwind(u, Tb, NJ, T0, C, tend)
% Returns the temperature distribution in a 1D system with initial
% condition T0 that is convected with velocity u. The temperature
% distribution is numerically calculated using first order upwind method
% with explicity Euler time discretization.

% Initalize temperature vectors
if length(T0) ~= NJ
    error("Length of initial temperature is not " + num2str(NJ));
end
T = T0;
Tnew = zeros(1,NJ);

% Time step
dt = C / (u*NJ);
t = 0;

% Calculate temperatue
while t <= tend
    for i = 1:1:NJ-1
        Tnew(i) = T(i)*(1+C) - C*T(i+1);
    end
    Tnew(NJ) = T(NJ)*(1+C) - C*Tb;
    T = Tnew;
    t = t + dt;
end
