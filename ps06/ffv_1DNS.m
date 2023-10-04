function f = ffv_1DNS(U, dx, gamma, Re, Pr)
% Viscous flux functions for nondimensionalized 1D Navier-Stokes
% equation
%
% (0  )

u = U(2,:)./U(1,:);
du = zeros(size(u));
du(2:end-1) = 1/(2*dx)*(u(3:end) - u(1:end-2));
du(1) = 1/(2*dx)*(u(1) - u(2));
du(end) = 1/(2*dx)*(u(end) - u(end-1));
p = (gamma-1) * (U(3,:) - 1/2*U(2,:)) .* u;
T = p./U(1,:);
dT = zeros(size(u));
dT(2:end-1) = 1/(2*dx)*(T(3:end) - T(1:end-2));
dT(2:end-1) = 1/(2*dx)*(T(3:end) - T(1:end-2));
dT(1) = 1/(2*dx)*(T(1) - T(2));
f = [zeros(1,length(du)); ...
     4/(3*Re).*du; ...
     1/Re*(2/3.*u.*du + gamma/(gamma-1)/Pr.*dT)];
end