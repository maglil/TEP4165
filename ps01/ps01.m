clc; close all; clear all

% Physical parameters
alpha = 1e-1; % Thermal diffusion coefficient

% System geometry
L1 = 1;
L2 = 6;

% Inital and boundary conditions
T1 = 263.15;
T2 = 393.15;
Ti = 293.15;

% Discretization
t = [0,5,10,15,20,25];
x = linspace(L1,L2,101);

f = figure; hold on

nstop = zeros(length(t)); % store truncation of each sum
for i = [1:1:length(t)]
    tval = t(i);
    [T,n(i)] = Texact(x,tval, alpha, T1, T2, Ti, L1, L2);
    plot(x,T, 'LineWidth',1)
end

% Figure properties
title('Temperature distrbution at different time points')
xlabel('Position (m)')
ylabel('Temperature (K)')
legendArray = strcat('t=',string(num2cell(t)));
h = legend(legendArray);
set(h, 'Location', 'NorthWest')
xlim([0.8,6.2]);
grid on

saveas(f,'ps01-1.png')
writematrix(n, "nstop.csv")
writematrix(t, "times.csv")
