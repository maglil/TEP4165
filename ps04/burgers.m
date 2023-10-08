function u = burgers(NJ, u0, nu, s, tend)

% order
n = 1; % padding for periodic boundary conditions

% Initalize vectors
u = zeros(1,NJ+2*n);
unew = u;

u(1+n:end-n) = u0; % Inital conditions

u(1:1+n) = u(NJ:NJ+n); % Periodic boundary conditions
u(NJ+n:end) = u(1+n:1+2*n);

% cellsize
dx = 1/NJ;

% timestep
dt = s / (max(abs(u0))/dx + 2*nu/dx^2);
t = 0; 

%For each timestep
while t < tend    
    fw = godunov_burgers(u(n), u(n+1)); % First western flux
    % For each cell
    for j = 1+n:NJ+n
        % Calculate flux
        fe = godunov_burgers(u(j), u(j+1)); % use upwind method
        
        % Solve discrete Burgers equation
        unew(j) = u(j) - dt/dx * (fe - fw) + nu*dt/dx^2*(u(j+1) + u(j-1) - 2*u(j));
        fw = fe; % Use eastern for next iteration step
    end
    unew(1:1+n) = unew(NJ:NJ+n); % Periodic boundary conditions
    unew(NJ+n:end) = unew(1+n:1+2*n);
    u = unew;
    t = t + dt;
end

u = u(1+n:NJ+n); % trim periodic padding
