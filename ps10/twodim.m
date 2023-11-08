% MATLAB program to load files temperature.dat, x.dat, y.dat and to plot T(x,y)
% TEP4165 Computational Heat and Fluid Flow, NTNU

clear all
close all
clc

%---Load numerical files
load 'temperature.dat'
load 'x.dat'
load 'y.dat'
load 'nusselt.dat'
load 'norms.dat'
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

figure(3)
plot(x(2:end-1),nusselt(2:end-1,1), x(2:end-1),nusselt(2:end-1,2), 'LineWidth', 1)
legend('Top wall', 'Bottom wall')
xlabel('x')
ylabel('Nusselt no.')
title('Nusselt number along wall')

figure(4)
semilogy(0:1:length(norms)-1, norms, 'LineWidth', 1)
xlabel('Iterations')
ylabel('2-norm')
title('Convergence behavior, u=0.05 m/s. Refined grid')