function f = flux_faces(F,method)
% Approximate fluxes at cell faces using cell center values
    if method == "central"
        f = 1/2*(F(:,1:end-1) + F(:,2:end));
    end
end