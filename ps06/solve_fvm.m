function U = solve_fvm(R, U0, ffc, ffv, tend, dx, alpha, gamma, Re, Pr, L, t_solver, approx)
% Solves the ODE dU/dt = R where the residual is determined using a finite
% volume approach.
% R Function that returns the residual from applying the FVM to a
% conservation law.
% U0 : Inital conditon
% ffc, ffv : flux function for convective and viscous fluxes 

% Notation:
% F (lowercase) for flux in cell centers
% f (uppercase) for flux on cell faces

dt = timestep(alpha, Pr, Re, gamma, U0, L, dx)
r = dt/dx;

U = U0;

if t_solver == "eeuler" % explicit Euler
    t = 0;
    while t < tend
        % Convective and viscous fluxes at cell centers
        Fc = ffc(U,gamma);
        Fv = ffv(U, dx, gamma, Re, Pr);
        % Approximate fluxes at cell faces
        fv = flux_faces(Fv, approx);
        fc = flux_faces(Fc, approx);
        
        U(:,2:end-1) = U(:,2:end-1) + r*R(fc, fv);
        U(:,end) = U(:,end-1);
        t = t+dt;
    end
end