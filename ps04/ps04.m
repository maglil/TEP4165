%% 2b
clc; clear; close all

% Physical parameters
nu = 0;

% Spatial discretization
NJ = 100;
dx = 1/NJ;
x = dx/2:dx:1-dx/2; % cell midpoint positions

% Inital condition
u0 = -1 -3*cos(2*pi*x);

% Stability criterion
Cmax = 1 ;
s = Cmax;

% Simulation length
%tend = 0.1;  run = 'test';
tend = 0.1:0.1:1; run = '';
%tend = 0.01:0.01:0.1; run = 'sf' % show shock formation

figure
hold on
for t = tend
    u = burgers(NJ,u0,nu,s,t);
    plot(x,u,'LineWidth',1.5)
end

% Graphics design
title("Inviscid Burgers equation, upwind method. " + run)
legendArray = strcat('t=',string(num2cell(tend)));
legend(legendArray)
box on
xlabel('x')
ylabel('u')

saveas(gcf,"ps04-2b" + run + ".png")

%% 1c-testing
clc; close all

% Physical parameters
nu = 0.02;

% Stability condition
s = 0.7;

% Simulation length
tend = 0.1;

u = burgers(NJ,u0,nu,s,t);

figure
plot(x,u,'LineWidth',1.5)

%% 1c-2
clc; close all

% Physical parameters
nu = 0.02;

% Stability condition
s = 0.7;

tend = 0.1:0.1:1;

figure
hold on

for t = tend
    u = burgers(NJ,u0,nu,s,t);
    plot(x,u,'LineWidth',1.5)
end

title(strcat('Burgers equation, upwind method. $\nu=$', num2str(nu)),'Interpreter','latex')
legendArray = strcat('t=',string(num2cell(tend)));
legend(legendArray)
xlabel('x')
ylabel('u')
box on

saveas(gcf,'ps04-2c.png')



