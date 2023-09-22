clc; clear; close all

nu = 0;

NJ = 100;
dx = 1/NJ;
x = [dx/2:dx:1-dx/2];
u0 = -1 -3*cos(2*pi*x);


Cmax = 1 ;
s = Cmax;

%tend = 1;
tend = 0.01:0.01:0.1;

figure
hold on
for t = tend
    t
    u = upwind(NJ,u0,nu,s,t);
    plot(x,u)
end