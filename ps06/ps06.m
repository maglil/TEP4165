% Compute the shock profile across a shock using NS.

% The state variable is U = [rho, rho*u, rho*E)

% Inital condition

NJ = 
U0 = %n x NJ
tend = 1

U = solve_fvm(@residual_1DNS, U0, @ffc_1DNS, @ffv, tend, dx, 