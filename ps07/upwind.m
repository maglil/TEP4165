function F = upwind(f,u, df, i)
% Calculate flux at cell faces using upwind method
% Face on right/upper face has same index as cell center
% Upper row, right column is not set
% f flux in i-direction in cell center
% u conserved variable
% df function handle to derivative of flux function 
    [NI,NJ] = size(f);
    F = zeros(NI,NJ);
    a = zeros(NI,NJ);
    

    charVel2 = @(fR,fL,uR,uL) charVel(fR,fL,uR,uL,df);

    if i == 1 % faces in x-direction      
        % Calculate characterisitic speed
        a(1:NI-1,:) = arrayfun(charVel2, f(2:NI,:), f(1:NI-1,:), u(2:NI,:), u(1:NI-1,:));
        % Calculate face flux using upwind method
        F(1:NI-1,:) = 1/2*(f(1:NI-1,:) + f(2:NI,:) - abs(a(1:NI-1,:)) .* (u(2:NI,:)-u(1:NI-1,:)) ) ;
    end

    if i == 2 % faces in y-direction      
        % Calculate characterisitic speed
        a(:,1:NI-1) = arrayfun(charVel2, f(:,2:NI), f(:,1:NI-1), u(:,2:NI), u(:,1:NI-1));
        % Calculate face flux using upwind method
        F(:,1:NI-1) = 1/2*(f(:,1:NI-1) + f(:,2:NI) - abs(a(:,1:NI-1)) .* (u(:,2:NI)-u(:,1:NI-1)) ) ;
    end
end