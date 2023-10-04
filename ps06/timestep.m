function dt = timestep(alpha, Pr, Re, gamma, U0, L, dx)

u = U0(2,1)/U0(1,1);
p = (gamma-1)*(U0(3,1)-1/2*U0(2,1)*u);
c = sqrt(gamma*p/U0(1,1));
beta = max(4/3, gamma/Pr);
mu_rho = Re*u*L; % mu/rho
dt = alpha*min(2*beta*mu_rho*(abs(u)+c)^2, dx^2/(2*beta*mu_rho));
end