function [T,n]= Texact(x, t, alpha, T1, T2, Ti, L1, L2)
    
	m = 200; %maximum iteration steps
		
    Ts = T1 + (T2-T1)/(L2-L1)*(x-L1); % steady state solution
    Tt = 0; % time dependent solution
    
    n = 1;
    while n<m
        temp =  2/(pi*n) * ((T2-Ti)*(-1)^n - (T1-Ti)) * sin(n*pi*(x-L1)/(L2-L1)) * exp(-alpha*n^2*pi^2/(L2-L1)^2*t) ;    
        Tt = Tt + temp;
        % Stop iteration if term in sum is less than machine precision 
        if max(abs(temp)) < eps 
            break
        end
        n = n+1;
    end
    T = Ts+Tt;
end