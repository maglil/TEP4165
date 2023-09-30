function U = euler1D(U0, tend, x1, x2, NJ, gamma, r, method)

dx = (x2-x1)/NJ;

dt = r*dx;
U = U0;
t = 0;

if method == "rusanov"
    while t<tend
        fe = rusanov(U,@flux_cons_gas,gamma); % find fluxe at cell faces       
        U = U - r*(fe(:,2:NJ+1) - fe(:,1:NJ)); % time step
        t = t+dt;
    end
end

rinv = 1/r;
if method == "lax-friedrich"
    while t<tend
        fe = lax_friedrich(U, @flux_cons_gas,gamma, rinv);
        U = U - r*(fe(:,2:NJ+1) - fe(:,1:NJ));
        t = t+dt;
    end
end
