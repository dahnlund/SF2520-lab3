% Computer exercise 3 Part 1, David Ahnlund Emil Gestsson
clc, clear;

N = 10; %Discretization resolution in x
M = 100000; %Discretization resolution in t
Lx = 1;
T = 2;
d = 0.35;
a = 1.2;
dx = Lx/N;
dt = T/M;

%% A
A = d*1/dx^2 * spdiags([1*ones(N-1,1) (-2*ones(N-1,1)) 1*ones(N-1,1)], -1:1, N-1, N-1);

%Adjust for Neumann boundary condition
A(end,end) = d*1/dx^2 * (-2/3);
A(end,end-1) = d*1/dx^2 * (2/3);

b = @(t) d/(dx^2)*[ sin(pi*t/a) * (t<=a) ;zeros(N-2,1)];

% Initial condition
u0 = zeros(N-1,1);

saved_u = zeros(N-1,M+1);
saved_u(:,1) = u0;

%Explicit Euler
uk = u0;
t = 0:dt:T;
for n = 2:length(t)
    u_new = uk + dt*(A*uk + b(t(n)));
    saved_u(:,n) = u_new;
    uk = u_new;
end

x = dx:dx:Lx-dx;
surf(t,x,saved_u)
shading interp

