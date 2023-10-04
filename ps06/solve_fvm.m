function U = solve_fvm(R, U0, ffc, ffv, tend, dx, s, method, approx)
% Solves the ODE dU/dt = R where the residual is determined using a finite
% volume approach.
% R Function that returns the residual from applying the FVM to a
% conservation law.

% Notation:
% f (lowercase) for flux in cell centers
% F (uppercase) for flux on cell faces

dt =
r = dt/dx;

U = U0;

if method == "eeuler" % explicit Euler
    t = 0;
    while t < tend
        % Convective and viscous fluxes at cell centers
        fc = ffc(U);
        fv = ffv(U);
        % Approximate fluxes at cell faces
        Fv = approx(fv);
        Fc = approx(fc);
        U = Un + r*R(Fc, Fv);
        t = t+dt;
    end
end