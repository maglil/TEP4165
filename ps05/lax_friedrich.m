function fe = lax_friedrich(U,ff,gamma,rinv)
% Computes the flux at cell faces using the Lax-Friedrich method
% and flux function ff and given conserved vector U'. r is dx/dt.

%gamma = 1.4;

% Initalize vector
% fe : flux across cell face 
% NJ+1 faces
[M,NJ] = size(U);
fe = zeros(M,NJ+1);

% Calculate value of flux function at cell centers
f = ff(U,gamma);

fe(:,2:NJ) = 1/2*( (f(:,1:NJ-1) + f(:,2:NJ)) - ...
    rinv.*(U(:,2:NJ) - U(:,1:NJ-1) ) );
fe(:,1) = f(:,1); 
fe(:,NJ+1) = f(:,NJ);
end