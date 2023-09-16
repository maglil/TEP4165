% TEP4165 Problem set 3
close all; clear all; clc;
part2h = 1;

% Physical parameters
u = -0.25;

% Discretization parameters
NJ = 100;
C = -0.9;
T0 = zeros(1,NJ);

% Simulation time
tend = 1.2;

% Inital conditions
if ~part2h
    T0(1:(NJ/2)-1) = 1000;
    T0(NJ/2:end) = 200;
end

% Inital conditions part 2h
if part2h
    dx = 1/NJ;
    x = [dx/2:dx:1];
    T0 = cos(2*pi*x);
end

% Boundary value
if ~part2h
    Tb = 200;
end
if part2h
    Tb = 1;
end

T = upwind(u, Tb, NJ, T0, C, tend);

dx = 1/NJ;
x = [dx/2:dx:1];
plot(x, T, 'LineWidth',1.5)

if ~part2h
    ylim([0,1100])
elseif part2h
    ylim([-1.1,1.1]);
end
grid on
xlabel("Position")
ylabel("Temperature")
title(["Temperature distribution at t=" + num2str(tend) + ".","Upwind method, explicit Euler. C = " + num2str(C)])
hold on

% Exact solution
if ~part2h
    j = uint8((0.5+u*tend)/dx + 0.5);
    Texact = ones(1,NJ)*1000;
    Texact(j:end) = 200;
end

% part 2h
if part2h 
    Texact = cos(2*pi*(x-u*tend));
    j = uint8((1+u*tend)/dx + 0.5);
    Texact(j:end) = 1;
end

plot(x,Texact,'--','LineWidth',1.5)
legend('Upwind','Exact','Location','East')
if ~part2h
    saveas(gcf,"T" + num2str(tend) + "C" + num2str(C) + ".png")
end
if part2h
    saveas(gcf,"cosT" + num2str(tend) + "C" + num2str(C) + ".png")
end