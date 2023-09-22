%% 2b

clc; clear; close all

nu = 0;

NJ = 100;
dx = 1/NJ;
x = [dx/2:dx:1-dx/2];
u0 = -1 -3*cos(2*pi*x);

Cmax = 1 ;
s = Cmax;

%tend = 1;
tend = 0.1:0.1:1;
%tend = 0.01:0.01:0.1;

figure
hold on
for t = tend
    u = upwind(NJ,u0,nu,s,t);
    plot(x,u,'LineWidth',1.5)
end
title('Inviscid Burgers equation, upwind method')
legendArray = strcat('t=',string(num2cell(tend)));
legend(legendArray)
box on
xlabel('x')
ylabel('u')
saveas(gcf,'ps04-2b.png')

%% 1c-1
clc; close all

nu = 0.02;
s = 0.7;
tend = 0.1;

figure
plot(x,u,'LineWidth',1.5)

%% 1c-2
clc; close all

s = 0.7;
tend = 0.1:0.1:1;

figure
hold on

for t = tend
    u = upwind(NJ,u0,nu,s,t);
    plot(x,u,'LineWidth',1.5)
end

title(strcat('Burgers equation, upwind method. $\nu=$', num2str(nu)),'Interpreter','latex')
legendArray = strcat('t=',string(num2cell(tend)));
legend(legendArray)
xlabel('x')
ylabel('u')
box on

saveas(gcf,'ps04-2c.png')

%%



