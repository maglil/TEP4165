function T = heateq1d(L1, L2, NJ, T0, r, tend, alpha, T1, T2)
% heateq1d Solves the 1D heat equation with given boundary and inital
% conditions
%

dx = (L2-L1)/NJ;
dt = r*dx^2/alpha;

%Initalize vectors
t = 0;
T = T0;
assert(length(T0) == NJ, "Length of inital temperature vector T0 does not match number of cells NJ");
Tnew  = zeros(1,NJ);

while t<tend
    % Update new endpoints
    Tnew(NJ) = T(NJ) + r*( -3*T(NJ)+T(NJ-1)+2*T2 );
    Tnew(1) = T(1) + r*( -3*T(1)+T(2)+2*T1 );
    % Update interior points
    for i = 2:1:NJ-1
        Tnew(i) = T(i) + r*( -2*T(i)+T(i-1)+T(i+1) );
    end
    % Assign new to old
    T = Tnew;
    % Update time
    t = t+dt;
end
