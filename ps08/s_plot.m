% MATLAB program to load file res.dat and to plot T(x); 
% TEP4165 Computational Heat and Fluid Flow, NTNU; 

load 'res.dat'; 
x=res(:,1);
y=res(:,2);
figure(1);
plot(x,y,'k-+'); 
xlabel('x [m]'); 
ylabel('T [K]'); 
npi = length(x);
title(sprintf('Steady temperature profile computed with %g grid points',npi));
