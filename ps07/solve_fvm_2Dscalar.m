function u = solve_fvm_2Dscalar(x,y,NI,NJ, sscond)
% x vector of xmin,xmax
% y vector of ymin,ymax
% NI, NJ number of grid points in x and y directions between xmin,xmax,
% ymin,ymax

% Cell size, grid spacing
dx = (x(2)-x(1))/NI;
dy = (y(2)-y(1))/NJ;


% F1, F2 fluxes at cell faces in x,y directions 
while du < sscriterion
    @f1(u) = u.^2/2;
    @f2(u) = u;
    F2 = f2(u)
    F1 = f1(u);
    R = (F1(2:NJ+1,:) - F1(1:NJ,:))/dy + (F2(:,2:NJ+1) - F2(:,1:NJ))/dx;
    u = u - dt*R; 


dx
dy

u = 0;

end