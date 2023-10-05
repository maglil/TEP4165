function R = residual_1DNS(Fc, Fv)
% Residual after applying FVM to 1D NS equation
% Fc convective flux
% Fv viscous flux

c_size = size(Fc);
v_size = size(Fv);

assert(all(c_size == v_size), 'size of fluxes do not match')
assert(length(c_size) == 2, 'Flux vector must be 2D')

R = ( Fv(:,2:end) - Fv(:,1:end-1)) - (Fc(:,2:end) - Fc(:,1:end-1));
end