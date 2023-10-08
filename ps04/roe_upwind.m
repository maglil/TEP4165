function f = roe_upwind(u1, u2)

f = 1/2 * (u1^2/2 + u2^2/2 - abs(1/2*(u1 + u2))*(u2 - u1) );
