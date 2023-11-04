% Computer exercise 3 Part 1, David Ahnlund Emil Gestsson
clc, clear; global u0_value d;

N = 100; %Discretization resolution in x

Lx = 1;
T = 2;
d = 0.35;
a = 1.2;

M = d*2*T*N^2/(Lx^2);  %Discretization resolution in t (according to theory)

%Equivalent stability condition: f(\Delta x, \Delta t) = \frac{\Delta t}{(\Delta x)^2} < \frac{1}{2d} = c

dx = Lx/N;
dt = T/M;

% dirichlet condition at u(0, t)
u0_value= @(t) sin(pi*t/a) .* (t<=a);

% add values at the boundaries
add_bounds = @(t, u) [u0_value(t); u; (4*u(end, :)-u(end-1, :)) / 3];

%% A
[u0, A, b, dudt] = create_system(N, dx);

saved_u = zeros(N-1,M+1);
saved_u(:,1) = u0;

%Explicit Euler
uk = u0;
t = 0:dt:T;
for n = 2:length(t)
    u_new = uk + dt*dudt(t(n), uk);
    saved_u(:,n) = u_new;
    uk = u_new;
end

saved_u = add_bounds(t, saved_u);

x = 0:dx:Lx;
surf(t,x,saved_u)
shading interp
zlabel("Temperatur")
xlabel("time, t")
ylabel("Rod position x")

%% Plot at specific tau

tau = 1.1;
u = saved_u(:,t==tau);

plot(x,u)
xlabel("x")
title("Plot of u at specific \tau")
legend("\tau = "+string(tau))

%% Plot at boundaries, as a function of time

u_leftend = saved_u(1, :);
u_rightend = saved_u(end, :);

plot(t,u_leftend)
hold on
plot(t, u_rightend)
xlabel("t")
ylabel("Temperature")
title("Plot of u at ends (through time)")
legend("x = 0", "x = " + string(Lx))

%% Using matlab functions:
%% ODE23

N_list = [100 200 400];

for i = 1:length(N_list)
    N = N_list(i);
    dx = Lx/N;
    
    [u0, A, b, dudt] = create_system(N, dx);
    options = odeset(RelTol=1e-4);

    tic;
    [t, u_ode23] = ode23(dudt, [0 2], u0, options);
    time = toc;
    u_ode23 = add_bounds(t', u_ode23');

    fprintf("ODE23: #Time steps = %.0d, for N = %.0d, CPU-time: %.08f seconds\n", length(t), N, time)
    
    if i == 1
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
    
    [u0, A, b, dudt] = create_system(N, dx);
    options = odeset(RelTol=1e-4);

    tic;
    [t, u_ode23s] = ode23s(dudt, [0 2], u0, options);
    time = toc;

    u_ode23s = add_bounds(t', u_ode23s');

    fprintf("ODE23s: #Time steps = %.0d, for N = %.0d, CPU-time: %.08f seconds\n", length(t), N, time)
    
    if i == 1
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
    
    [u0, A, b, dudt] = create_system(N, dx);
    options = odeset("Jacobian",A,RelTol=1e-4);
    
    tic;
    [t, u_ode23sJ] = ode23s(dudt, [0 2], u0, options);
    time = toc;
    
    u_ode23sJ = add_bounds(t', u_ode23sJ');

    fprintf("ODE23sJ: #Time steps = %.0d, for N = %.0d, CPU-time: %.08f seconds\n", length(t), N, time)
    
    if i == 1
        x = 0:Lx/N:Lx;
        surf(t,x,u_ode23sJ)
        shading interp
    end
end
fprintf("\n\n")

%% Avoid repeat code
function [u0, A, b, dudt] = create_system(N, dx)
    global u0_value d;

    u0 = zeros(N-1,1);
    A = d*1/dx^2 * spdiags([1*ones(N-1,1) (-2*ones(N-1,1)) 1*ones(N-1,1)], -1:1, N-1, N-1);

    %Adjust for Neumann boundary condition
    A(end,end) = d*1/dx^2 * (-2/3);
    A(end,end-1) = d*1/dx^2 * (2/3);
    
    b = @(t) d/(dx^2)*[ u0_value(t); zeros(N-2,1)];
    dudt = @(t,u) A*u+b(t);
end