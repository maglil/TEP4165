%% TEP 4280, Problem set 5

%set(groot, 'DefaultLineWidth', 1.5);

%% Task 1c
clear; clc; close all

gamma = 1.4;

NJ = 800;
U0 = zeros(3,NJ);

U0(:,1:NJ/2) = repmat([1;0;2500],1,NJ/2);
U0(:,NJ/2+1:NJ) = repmat([1;0;0.025],1,NJ/2);

tend = 0.024;
x1 = -1;
x2 = 1;

r = 0.0187;

U = euler1D(U0, tend, x1, x2, NJ, gamma, r, "rusanov");

% Plot
dx = (x2-x1)/NJ;
x = -1:dx:1-dx;

rho = U(1,:);
u = U(2,:)./rho;
p = U(3,:) - 1/2*rho.*u.^2;

key = {"LineWidth", "Color"};
val = {1.5,'k'};

h1=plot(x,rho);
set(h1,key,val)
xlabel("\rho")
ylabel("x")
title("Density")
grid on
saveas(gcf,'rho.png')

figure
h2 = plot(x,u);
set(h2,key,val)
xlabel("u")
ylabel("x")
title("Velocity")
grid on
saveas(gcf, 'u.png')

figure
h3 = plot(x,p);
set(h3,key,val)
xlabel("p")
ylabel("x")
title("Pressure")
grid on
saveas(gcf, 'p.png')

%% Task 1d
close all; clear; clc

gamma = 1.4;

NJ = 800;
U0 = zeros(3,NJ);

U0(:,1:NJ/2) = repmat([1;0;2500],1,NJ/2);
U0(:,NJ/2+1:NJ) = repmat([1;0;0.025],1,NJ/2);

tend = 0.024;
x1 = -1;
x2 = 1;


r1 = 0.0187;
r2 = 0.0001;
tic
U1 = euler1D(U0, tend, x1, x2, NJ, gamma, r1, "rusanov");
toc
tic
U2 = euler1D(U0, tend, x1, x2, NJ, gamma, r2, "rusanov");
toc

% Plot
dx = (x2-x1)/NJ;
x = -1:dx:1-dx;

rho1 = U1(1,:);
u1 = U1(2,:)./rho1;
p1 = U1(3,:) - 1/2*rho1.*u1.^2;

rho2 = U2(1,:);
u2 = U2(2,:)./rho2;
p2 = U2(3,:) - 1/2*rho2.*u2.^2;

key = {"LineWidth"};
val = {1.5};

h1=plot(x,rho1, x, rho2);
set(h1,key,val)
xlabel("\rho")
ylabel("x")
legend("\rho_1", "\rho_2")
title("Density")
grid on
saveas(gcf,'rho_dt.png')

figure
h2 = plot(x,u1, x, u2);
set(h2,key,val)
xlabel("u")
ylabel("x")
legend("u_1", "u_2")
title("Velocity")
grid on
saveas(gcf, 'u_dt.png')

figure
h3 = plot(x,p1, x, p2);
set(h3,key,val)
xlabel("p")
ylabel("x")
legend("p_1", "p_2")
title("Pressure")
grid on
saveas(gcf, 'p_dt.png')
