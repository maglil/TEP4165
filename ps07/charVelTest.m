% Test script for charvel
% runtests('charVelTest')
clc

%% Test 1 : uR = uL, df(u)=1 returned
df = @(u) 1;
fR = 1;
fL = 2;
uR = 1;
uL = 1;
a = charVel(fR, fL, uR, uL, df);
assert(a == 1, 'Failed uR=uL, df=const')

%% Test 2 : uR = uL, df(u)=u returned
df = @(u) u;
fR = 1;
fL = 2;
uR = 2;
uL = 2;
a = charVel(fR, fL, uR, uL, df);
assert(a == 2, 'Failed uR=uL, df=u')

%% Test 3 : uR != uL, (fR-fL)/(uR-uL) returned
df = @(u) u;
fR = 4;
fL = 2;
uR = 2;
uL = 1;
a = charVel(fR, fL, uR, uL, df);
assert(a == 2, 'Failed test 3')