function f = ffv_1DNS(U,gamma, Re, Pr)
% Viscous flux functions for nondimensionalized 1D Navier-Stokes
% equation
%
% (0  )

u = U(2,:)./U(1,:);
p = (gamma-1) * (U(3,:) - 1/2*U(2,:)) * u;
T = p./U(1,:);
f = [0; ...
     4/(3*Re).*u; ...
     1/Re*(2/3.*u.^2 + gamma/(gamma-1)/Pr.*T)];
end