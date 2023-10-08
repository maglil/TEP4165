function f = godunov_burgers(u1,u2)

if u1 == u2
    f = u1^2/2;
end
if u1<u2
    if sign(u1*u2) <= 0
        f = 0;
    elseif u1>0
        f = u1^2/2;
    else
        f = u2^2/2;
    end
end
if u1>u2
    if u1>abs(u2)
        f = u1^2/2;
    else
        f = u2^2/2;
    end
end