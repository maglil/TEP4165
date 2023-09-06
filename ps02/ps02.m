clear all, clc, close all
addpath("../ps01")

% Set inital conditions

NJ = 50;
alpha = 1e-2;
L1 = -3;
L2 = 5;
Ti = 305.15;
T0 = Ti*ones(1,NJ);
tend = 40;
T1 = 253.15;
T2 = 443.15;
r = 0.25;

T = heateq1d(L1, L2, NJ, T0, r, tend, alpha, T1, T2);

% Create arrays including boundary values and x position at cell centers
dx = (L2-L1)/NJ;
x = zeros(1,NJ+2);
x(2:end-1) = linspace(L1+dx,L2-dx, NJ);
x(1) = L1;
x(end) = L2;
Tbv = zeros(1,NJ+2);
Tbv(2:end-1) = T;
Tbv(1) = T1;
Tbv(end) = T2;

plot(x,Tbv);
x = L1 + (L2-L1)/2+dx;
[Tx, n] = Texact(x, tend, alpha, T1, T2, Ti, L1, L2);


