% Computer exercise 3 Part 1, David Ahnlund Emil Gestsson
clc, clear;

N = 100; %Discretization resolution in x

Lx = 1;
T = 2;
d = 0.35;
a = 1.2;

M = d*2*T*N^2/(Lx^2);  %Discretization resolution in t (according to theory)

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

%Apply boundaries to solution
uN = 1/3*(4*saved_u(end,:)-saved_u(end-1,:));
u_initial = sin(pi*t/a) .* (t<=a);
saved_u = [u_initial;saved_u;uN];

x = 0:dx:Lx;
surf(t,x,saved_u)
shading interp

%% Plot at specific tau

tau = 1.1;
u = saved_u(:,t==tau);

plot(x,u)
xlabel("x")
title("Plot of u at specific \tau")
legend("\tau = "+string(tau))

%% Using matlab functions:
%% ODE23

N_list = [100 200 400];

for i = 1:length(N_list)
    N = N_list(i);
    dx = Lx/N;
    
    u0 = zeros(N-1,1);

    A = d*1/dx^2 * spdiags([1*ones(N-1,1) (-2*ones(N-1,1)) 1*ones(N-1,1)], -1:1, N-1, N-1);

    %Adjust for Neumann boundary condition
    A(end,end) = d*1/dx^2 * (-2/3);
    A(end,end-1) = d*1/dx^2 * (2/3);
    
    b = @(t) d/(dx^2)*[ sin(pi*t/a) * (t<=a) ;zeros(N-2,1)];

    dudt = @(t,u) A*u+b(t);
    options = odeset(RelTol=1e-4);

    tic;
    [t, u_ode23] = ode23(dudt, [0 2], u0, options);
    time = toc;
    u_ode23 = u_ode23';
    
    %Apply boundaries to solution
    uN = 1/3*(4*u_ode23(end,:)-u_ode23(end-1,:));
    u_initial = (sin(pi*t/a) .* (t<=a))';
    u_ode23 = [u_initial;u_ode23;uN];
    fprintf("ODE23: #Time steps = %.0d, for N = %.0d, CPU-time: %.08f seconds\n", length(t), N, time)
    
    if N == 100
        x = 0:Lx/N:Lx;
        surf(t,x,u_ode23)
        shading interp
    end
end
fprintf("\n\n")
%% ODE23s

N_list = [100 200 400];

for i = 1:length(N_list)
    N = N_list(i);
    dx = Lx/N;
    
    u0 = zeros(N-1,1);

    A = d*1/dx^2 * spdiags([1*ones(N-1,1) (-2*ones(N-1,1)) 1*ones(N-1,1)], -1:1, N-1, N-1);

    %Adjust for Neumann boundary condition
    A(end,end) = d*1/dx^2 * (-2/3);
    A(end,end-1) = d*1/dx^2 * (2/3);
    
    b = @(t) d/(dx^2)*[ sin(pi*t/a) * (t<=a) ;zeros(N-2,1)];

    dudt = @(t,u) A*u+b(t);
    options = odeset(RelTol=1e-4);

    tic;
    [t, u_ode23s] = ode23s(dudt, [0 2], u0, options);
    time = toc;
    u_ode23s = u_ode23s';
    
    %Apply boundaries to solution
    uN = 1/3*(4*u_ode23s(end,:)-u_ode23s(end-1,:));
    u_initial = (sin(pi*t/a) .* (t<=a))';
    u_ode23s = [u_initial;u_ode23s;uN];
    fprintf("ODE23s: #Time steps = %.0d, for N = %.0d, CPU-time: %.08f seconds\n", length(t), N, time)
    
    if N == 100
        x = 0:Lx/N:Lx;
        surf(t,x,u_ode23s)
        shading interp
    end
end
fprintf("\n\n")
%% ODE23sJ

N_list = [100 200 400];

for i = 1:length(N_list)
    N = N_list(i);
    dx = Lx/N;
    
    u0 = zeros(N-1,1);

    A = d*1/dx^2 * spdiags([1*ones(N-1,1) (-2*ones(N-1,1)) 1*ones(N-1,1)], -1:1, N-1, N-1);

    %Adjust for Neumann boundary condition
    A(end,end) = d*1/dx^2 * (-2/3);
    A(end,end-1) = d*1/dx^2 * (2/3);
    
    b = @(t) d/(dx^2)*[ sin(pi*t/a) * (t<=a) ;zeros(N-2,1)];

    dudt = @(t,u) A*u+b(t);
    options = odeset("Jacobian",A,RelTol=1e-4);
    
    tic;
    [t, u_ode23sJ] = ode23s(dudt, [0 2], u0, options);
    time = toc;
    u_ode23sJ = u_ode23sJ';
    
    %Apply boundaries to solution
    uN = 1/3*(4*u_ode23sJ(end,:)-u_ode23sJ(end-1,:));
    u_initial = (sin(pi*t/a) .* (t<=a))';
    u_ode23sJ = [u_initial;u_ode23sJ;uN];
    fprintf("ODE23sJ: #Time steps = %.0d, for N = %.0d, CPU-time: %.08f seconds\n", length(t), N, time)
    
    if N == 100
        x = 0:Lx/N:Lx;
        surf(t,x,u_ode23sJ)
        shading interp
    end
end
fprintf("\n\n")