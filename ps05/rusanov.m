function fe = rusanov(U,ff, gamma)
% Function calculates the fluxes at cell edges using Rusanovs method (local
% Lax-Friedrich). 
%
% ff is flux function
% U = [rho, rho*u, rho*E)
%
% Boundary conditions implemented as constant u across cell face.
%

%gamma = 1.4;

% fe : flux across cell face 
% NJ+1 faces
[M,NJ] = size(U);
fe = zeros(M,NJ+1);

f = ff(U,gamma);

% Local speed of sound in each cell
p = (gamma-1)*(U(3,:) - 1/2 * U(2,:).^2./U(1,:));
c = ( gamma*p./U(1,:) ).^(1/2);

% local 'convection velocity' on each cell face
a = zeros(1,NJ+1);
u = U(2,:)./U(1,:);
a(2:NJ) = max(abs(u(1:NJ-1)) + c(1:NJ-1), abs(u(2:NJ)) + c(2:NJ)); % interior points
a(1) = abs(u(1)) + c(1); % boundary 
a(NJ) = abs(u(NJ)) + c(NJ);

% Calculate flux at cell edge using Rusanovs method
% Interior point
fe(:,2:NJ) = 1/2*(f(:,1:NJ-1) + f(:,2:NJ)) - ...
    1/2*a(2:NJ).*(U(:,2:NJ) - U(:,1:NJ-1) );
% boundaries
fe(:,1) = f(:,1); 
fe(:,NJ+1) = f(:,NJ);
