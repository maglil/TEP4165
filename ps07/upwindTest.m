% Test script for upwind
% runtests('upwindTest')
clc

u = [1 2 3; 4 5 6; 7 8 9 ];

f = u;
df = @(u) 1;
F = upwind(f,u,df, 1);

%% Test 1 
assert(F(1,1)==1);

%% Test 2
assert(F(1,2)==2);

%% Test 2 : last row not set
assert(F(3,1)==0);

%% Test 4 : 
f = u.^2; % characteristic velocity always positive f is uL^2
df = @(u) 2.*u;
F = upwind(f,u,df, 1);
assert(F(1,2) == 4) 

%% Test 5 : 
f = u.^2; % characteristic velocity always positive f is uL^2
df = @(u) 2.*u;
F = upwind(f,u,df, 1);
assert(F(2,1) == 16) 

%% Test 6
u = [4 5 6; 1 2 3; 7 8 9 ];
f = u;
df = @(u) 1;
F = upwind(f,u,df, 1);
assert(F(1,1)==4)

%% Test 7
u = [4 5 6; 1 2 3; 7 8 9 ];
f = u;
df = @(u) 1;
F = upwind(f,u,df, 1);
assert(F(1,2)==5)

%% Test 8
u = [1 2 3; -2 -3 -4; 7 8 9 ];
f = u.^2;
df = @(u) 2.*u;
F = upwind(f,u,df, 1);
F(1,2)
assert(F(1,2)==9)

%% Test 9
u = [1 2 3; -2 -3 -4; 7 8 9 ];
f = u.^2;
df = @(u) 2.*u;
F = upwind(f,u,df, 1);
F(1,3)
assert(F(1,3)==16)

%% Test 10
u = [1 2 3; 4 5 6; 7 8 9 ];
f = u.^2;
df = @(u) 2.*u;
F = upwind(f,u,df, 2);
F(1,3)
assert(F(2,2)==25)


