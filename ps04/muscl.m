function u = muscl(NJ, u0, nu, s, tend)

unew = zeros(1,NJ+4); % Include periodic boundaries
u(3:NJ+2) = u0;
u(1) = u0(NJ-1);
u(2) = u0(NJ);
u(NJ+3) = u0(1);
u(NJ+4) = u0(2);

dx = 1/NJ;
dt = s / (max(abs(u0))/dx + 2*nu/dx^2);
t = 0;

% minmod function
minmod = @(a,b) sign(a) * max(0, min( abs(a) , sign(a)*b ));

while t<tend
    for j = 3:NJ+2
        % Calculate fluxes
        uL(j) = u(j) + 1/2*minmod(u(j)-u(j-1),u(j+1)-u(j));
        uR(j) = u(j+1) - 1/2*minmod(u(j+2)-u(j+1),u(j+1)-u(j));
        un(j) = uR(j)^2/2 - uL(j)^2/2;
        % Riemann problem       
        unew(j) = u(j) - dt/dx*(uR(j)^2/2 - uL(j)^2/2);
    end
     % Update periodic boundary conditions
    unew(1) = unew(NJ+1);
    unew(2) = unew(NJ+2);
    unew(NJ+3) = unew(3);
    unew(NJ+4) = unew(4);

    u = unew;
    t = t + dt;
end
u = u(3:NJ+2);
