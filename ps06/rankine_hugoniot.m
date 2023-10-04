function [rho2, u2, p2] = rankine_hugoniot(rho1, u1, p1, gamma, Ma)

    f = 2/(gamma-1)*(1-1/Ma^2);
    rho2 = rho1*(1-f)^(-1);
    u2 = u1*(1-f)^(-1);
    p2 = p1*(1 + 2*gamma/(gamma+1)*(Ma^2-1));

end