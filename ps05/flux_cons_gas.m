function ff = flux_cons_gas(U, gamma)
    % flux function for conservative variables and ideal gas with
    % given ratio of heat capacaties.
    
    p = (gamma-1) * (U(3,:) - 1/2*U(2,:).^2./U(1,:) ) ;
    assert(min(p>=0)==1,"Negative pressure detected");
    
    ff1 = U(2,:);
    ff2 = U(2,:).^2./U(1,:) + p;
    ff3 = (U(3,:) + p).*U(2,:)./U(1,:);
    ff = [ff1;ff2;ff3];