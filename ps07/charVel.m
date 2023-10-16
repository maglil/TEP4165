function a = charVel(fR, fL, uR, uL, df)
% Computes the characteristic velocity

if uR ~= uL
    a = (fR-fL)/(uR-uL);
else
    a = df(uR);
end