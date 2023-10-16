function [u,X,Y,du] = solve_fvm_2Dscalar(xlim,ylim,NI,NJ, dt, sscond)
% x vector of xmin,xmax
% y vector of ymin,ymax
% NI, NJ number of grid points in x and y directions between xmin,xmax,
% sscond steady state condition
% dt timestep

%dt = 0.01;
%dt = 0.005; %N = 50
%dt = 0.002; %N = 100

% Cell size, grid spacing
dx = (xlim(2)-xlim(1))/NI;
dy = (ylim(2)-ylim(1))/NJ;
x = linspace(xlim(1), xlim(2), NI);
y = linspace(ylim(1), ylim(2), NJ);
[Y,X] = meshgrid(x,y); % x is the first dimension, that is rows, opposite to meshgrid convention
u = 0.75-2.*X;
% u = zeros(NI,NJ);

% flux functions
ff1 = @(u) u.^2/2;
df1 = @(u) u;
ff2 = @(u) u;
df2 = @(u) 1;
%assert(isa(df1,'function_handle'));

% Set boundary conditions
u(1,:) = 0.75;
u(NI,:) = -1.25;
u(:,1) = 0.75 - 2.*x;
un = u;

cnt = 0;
du_norm = 1;
du_norm0 = 1;
while du_norm > sscond*du_norm0
    cnt = cnt + 1;
    % Flux functions
    f1 = ff1(u);
    f2 = ff2(u);

    % F1, F2 fluxes at cell faces in x,y directions, respectively.
    F1 = upwind(f1, u, df1, 1);
    F2 = upwind(f2, u, df2, 2);

    F2(:,NJ) = f2(:,NJ); % Assume characteristic speed upwind on top boundary
    
    % Residual
    R = (F1(2:NI-1,2:NJ) - F1(1:NI-2,2:NJ))/dy + (F2(2:NI-1,2:NJ) - F2(2:NI-1,1:NJ-1))/dx;
    
    % time step
    un(2:NI-1, 2:NJ) = u(2:NI-1, 2:NJ) - dt*R; 
    
    du_norm = (dx*dy*sum(sum((un - u).^2)))^(1/2);
    
    if cnt==1
        du_norm0 = du_norm;
    end

    u = un;
    
    du(cnt) = du_norm;

    if du_norm <= sscond*du_norm0 % plot some params before exit
        du_norm
        cnt
    end

    % Maximum iterations
    if cnt>500
        "Maximum number of iterations reached."
        cnt
        du_norm
        break
    end
end