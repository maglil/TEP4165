clear all, clc, close all
addpath("../ps01")

%% Properties
% Material parameters
alpha = 1e-2;

% Discretization parameters
NJ = 50;
r = 0.25;
tend = 40;

% Geometry
L1 = -3;
L2 = 5;

% Inital conditions
Ti = 305.15;
T0 = Ti*ones(1,NJ);

% Boundary conditions
T1 = 253.15;
T2 = 443.15;


%% Task 2c and d

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

plot(x,Tbv, 'LineWidth', 1.5);
title("Temperature distribution for NJ = " + num2str(NJ) + ", tend = " + num2str(tend))
xlabel("Position")
ylabel("Temperature")
saveas(gcf,"t_" + num2str(tend) + ".png")

%% Task 2e
r = 0.6;
T = heateq1d(L1, L2, NJ, T0, r, tend, alpha, T1, T2);
Tbv(2:end-1) = T;

plot(x,Tbv, 'LineWidth', 1.5);
title("Stability criterion violoated, r=0.6")
xlabel("Position")
ylabel("Temperature")
saveas(gcf,"unstable" + ".png")

%% Task 2f and 2g

r = 0.25;
NJ = [50, 100, 200, 400];
error2 = zeros(1,length(NJ));
for i = 1:1:length(NJ)              
     % Set up arrays including boundary conditions
    dx = (L2-L1)/NJ(i);
    x = zeros(1,NJ(i)+2);
    x(2:end-1) = linspace(L1+dx,L2-dx, NJ(i));
    x(1) = L1;
    x(end) = L2;
    Tbv = zeros(1,NJ(i)+2);
    Tbv(1) = T1;
    Tbv(end) = T2;    
    % exact solution
    [Tx, n] = Texact(x, tend, alpha, T1, T2, Ti, L1, L2);
    % Numerical solution
    T0 = Ti*ones(1,NJ(i));
    T = heateq1d(L1, L2, NJ(i), T0, r, tend, alpha, T1, T2);
    Tbv(2:end-1) = T;
    
    error2(i) = 1/NJ(i)*(sum(T-Tx(2:end-1)))^(1/2);
    
    figure()
    plot(x,Tbv, 'LineWidth', 1.5)
    hold on 
    plot(x,Tx, 'LineWidth', 1.5)
    title("Numerical and anlytic solution for NJ = " + num2str(NJ(i)))
    xlabel("Position")
    ylabel("Temperature")
    legend("Exact", "Numerical")
    saveas(gcf,"numerical_analytic" + num2str(NJ(i)) + ".png") 
end

%% Task 1h
p = log(error2(1:3)./error2(2:4))/log(2);

