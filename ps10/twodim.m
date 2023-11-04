% MATLAB program to load files temperature.dat, x.dat, y.dat and to plot T(x,y)
% TEP4165 Computational Heat and Fluid Flow, NTNU

clear all
close all
clc

%---Load numerical files
load 'temperature.dat'
load 'x.dat'
load 'y.dat'
%---Plot the temperature
T=temperature';
[X,Y] = meshgrid(x,y);

figure(1)
surface(X,Y,T)
view(2)
daspect([1 1 1])
xlabel('x')
ylabel('y')
ylim([0,max(y)])
title('Temperature')
shading interp
lighting phong
colorbar
 
figure(2)
surfc(X,Y,T)
xlabel('x')
ylabel('y')
zlabel('T')
