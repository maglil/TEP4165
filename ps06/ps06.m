% Compute the internal profile across a shock using NS.
% The state variable is U = [rho, rho*u, rho*E), nondimensionalized.

%% Set up problem
clc, clear all, close all

% Physical parameters
Re = 1;
Pr = 0.75;
gamma = 1.4;

% Geometry
x1 = -25;
x2 = 25;
L = x2-x1;

% Space discretization
NJ = 100;
dx = NJ/L;

% Inital condition
Ma = 3;
U0 = zeros(3,NJ);
Ul = [1; sqrt(gamma)*Ma; 1/(gamma-1) + 1/2*gamma*Ma^2];
[rho2, u2, p2] = rankine_hugoniot(1, Ul(2), 1,gamma, Ma);
Ur = [rho2; rho2*u2; p2/(gamma-1) + 1/2*rho2*u2^2];
U0(:,1:NJ/2) = repmat(Ul,1,NJ/2);
U0(:,NJ/2+1:NJ) = repmat(Ur,1,NJ/2);

% Stability buffer
alpha = 0.125;

% Simulation duration
tend = 0.3;

%% Solve
clc
%L = 1? reference length
U = solve_fvm(@residual_1DNS, U0, @ffc_1DNS, @ffv_1DNS, tend, dx, alpha, gamma, Re, Pr, L, "eeuler", "central");
x = linspace(-25+dx,25-dx,100);
plot(x,U(1,:),x,U(2,:),x,U(3,:))
