function u = upwind(NJ, u0, nu, s, tend)

% Initalize vectors
unew = zeros(1,NJ+2); % Include periodic boundaries
u(2:NJ+1) = u0;
u(1) = u0(NJ);
u(NJ+2) = u0(1);

dx = 1/NJ;
dt = s / (max(abs(u0))/dx + 2*nu/dx^2);
t = 0;
while t<tend
    for j = 2:NJ+1
        fe = u(j)^2/2 + u(j+1)^2/2 - abs(1/2*(u(j) + u(j+1)))*(u(j+1) - u(j));
        fw = u(j-1)^2/2 + u(j)^2/2 - abs(1/2*(u(j-1) + u(j)))*(u(j) - u(j-1));
        unew(j) = u(j) - dt/(2*dx) * (fe-fw) + nu*dt/dx^2*(u(j+1) + u(j-1) - 2*u(j));        
        unew(1) = u(NJ+1);
        unew(NJ+2) = u(2); 
    end
    u = unew;
    t = t + dt;
end
u = u(2:NJ+1);
