% TEP4165 Problem set 3

% Physical parameters
u = -0.25;

% Discretization parameters
NJ = 100;
C = -0.75;
T0 = zeros(1,NJ);

% Simulation time
tend = 0.05;

% Inital conditions
T0(1:(NJ/2)-1) = 1000;
T0(NJ/2:end) = 200;

% Boundary value
Tb = 200;

T = upwind(u, Tb, NJ, T0, C, tend);

dx = 1/NJ;
x = linspace(0+dx/2,1-dx/2,NJ);
plot(x, T)