clc; clear all; close all;

NI = 100;
NJ = 100;
xlim = [0,1];
ylim = [0,1];

N = [25,50,100];
dts = 1./N/3; % time steps from stability condition.
% dts = [0.01,0.005, 0.004]; % set manually

for i = 1:length(N)
    NI = N(i);
    NJ = NI;
    dt = dts(i);
    [u,X,Y,du] = solve_fvm_2Dscalar(xlim,ylim,NI,NJ,dt,1e-4);
    contour(X,Y,u)
    title("N = " + num2str(NI))
    xlabel("x")
    ylabel("y")
    saveas(gcf,"ps07" + num2str(NI) + ".png")
    d{i} = du;
end
