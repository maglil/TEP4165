

% minmod function
minmod = @(a,b) sign(a) * max(0,min( abs(a) , sign(a)*b ));

% Calculate fluxes
uL = u(j) + 1/2*minmod(u(j)-u(j-1),u(j+1)-u(j));
uR = u(j+1) - 1/2*minmod(u(j+2)-u(j+1),u(j+1)-u(j));

% Riemann problem

unew(j) = u(j) - dt/dx(uR^2/2 - uL^2/2)