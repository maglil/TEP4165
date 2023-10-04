function f = ffc_1DNS(U,gamma)
% Convective flux functions for nondimensionalized 1D Navier-Stokes
% equation
%
% (rho*u; rho*u^2 + p; (rho*E + p)*u )

u = U(2,:)./U(1,:);
p = (gamma-1) * (U(3,:) - 1/2*U(2,:)) * u;
f = [U(2,:); ...
     U(1,:).*u + p; ...
     (U(3,:) + p).*u];
end