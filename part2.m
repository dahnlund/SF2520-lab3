% Computer exercise 3 Part 2, David Ahnlund Emil Gestsson

clear, clc;
%Coefficients
Lx = 12; Ly = 5; T_ext = 25;
%% A
N = 240;
h = Lx/N;
M = Ly/h;

F_func = @(x,y) 100*exp(-1/2 * (x-4).^2 - 4*(y-1).^2);

x = h:h:Lx-h;
y = h:h:Ly-h;
F = F_func(x',y);
Sx = 1/h^2*spdiags([1*ones(N-1,1) (-2*ones(N-1,1)) 1*ones(N-1,1)], -1:1, N-1, N-1);
Sy = 1/h^2*spdiags([1*ones(M-1,1) (-2*ones(M-1,1)) 1*ones(M-1,1)], -1:1, M-1, M-1);

%Boundary condition for x
Sx(1,1) = -2/(3*h^2); Sx(1,2) = 2/(3*h^2);
Sx(end,end) = -2/(3*h^2); Sx(end, end-1) = 2/(3*h^2);

%Boundary condition for y
Sy(end,end) = -2/(3*h^2); Sy(end, end-1) = 2/(3*h^2);

A = kron(speye(size(Sy)),Sx) + kron(Sy, speye(size(Sx)));

F(:,1) = F(:,1) + T_ext/h^2;

f = reshape(F, (N-1)*(M-1),1);

dt = 0.1;
u0 = T_ext * ones(size(f));

%Crank Nicolson
tau = 0:dt:40;
saved_u = zeros(length(u0),length(tau));
saved_u(:,1) = u0;

uk = u0;
for i = 2:length(tau)
    u_new = (speye(size(A)) - 1/2 * dt*A)\((speye(size(A)) + 1/2 * dt*A)*uk + dt*f);
    %u_new = (speye(size(A)) - 1 * dt*A)\(uk + dt*f);
    saved_u(:,i) = u_new;
    uk = u_new;
end

%Plot the solution
for frame = 1:1:40/dt
    tf = frame;

    y = 0:h:Ly;
    x = 0:h:Lx;

    u = reshape(saved_u(:,tf), (N-1), (M-1));
    u_y0 = T_ext * ones(N-1,1);
    u_M = 1/3*(4*u(:,end)-u(:,end-1));
    u = [u_y0 u u_M];   % Adding boundaries along y
    u_N = 1/3*(4*u(end,:)-u(end-1,:));
    u_x0 = 1/3*(4*u(1,:)-u(2,:));
    u = [u_x0;u;u_N];  %Adding boundaries along x

    mesh(y,x,u); % Update the plot
    view([-71.1 21.5877372262774]);
    zlim([25 70])
    xlabel("y")
    ylabel("x")
    
    drawnow; % Refresh the plot
end

%% C

u = saved_u(:,tau==40);
u = reshape(u, (N-1), (M-1));
u_y0 = T_ext * ones(N-1,1);
u_M = 1/3*(4*u(:,end)-u(:,end-1));
u = [u_y0 u u_M];   % Adding boundaries along y
u_N = 1/3*(4*u(end,:)-u(end-1,:));
u_x0 = 1/3*(4*u(1,:)-u(2,:));
u = [u_x0;u;u_N];  %Adding boundaries along x

fprintf("u(6,2,40) = %.04f\n\n", u(x==6, y==2))
   