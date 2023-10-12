% Computer exercise 3 Part 2, David Ahnlund Emil Gestsson

clear, clc;
%Coefficients
Lx = 12; Ly = 5; T_ext = 25;
%% A
N = 60;
h = Lx/N;
M = Ly/h;

F_func = @(x,y) 100*exp(-1/2 * (x-4).^2 - 4*(y-1).^2);

x = h:h:Lx-h;
y = h:h:Ly-h;
F = F_func(x',y);
Sx = -1/h^2*spdiags([1*ones(N-1,1) (-2*ones(N-1,1)) 1*ones(N-1,1)], -1:1, N-1, N-1);
Sy = -1/h^2*spdiags([1*ones(M-1,1) (-2*ones(M-1,1)) 1*ones(M-1,1)], -1:1, M-1, M-1);

%Boundary condition for x
Sx(1,1) = 2/(3*h^2); Sx(1,2) = -2/(3*h^2);
Sx(end,end) = 2/(3*h^2); Sx(end, end-1) = -2/(3*h^2);

%Boundary condition for y
Sy(end,end) = 2/(3*h^2); Sy(end, end-1) = -2/(3*h^2);

A = kron(speye(size(Sy)),Sx) + kron(Sy, speye(size(Sx)));

F(:,1) = F(:,1) + T_ext/h^2;

f = reshape(F, (N-1)*(M-1),1);

dt = 0.1;
u0 = T_ext * ones(size(f));

%Crank Nicholson
tau = 0:dt:40;
saved_u = zeros(length(u0),length(tau));
saved_u(:,1) = u0;

uk = u0;
for i = 2:length(tau)
   
    u_new = (eye + 1/2 * dt*A)\((eye + 1/2 * dt*A)*uk + dt*f);
    saved_u(:,i) = u_new;
    uk = u_new;

end

y = h:h:Ly-h;
x = h:h:Lx-h;
sol1 = reshape(saved_u(:,tau==40), (N-1), (M-1));
mesh(y,x,sol1)
xlabel("y")
ylabel("x")

%%
T = reshape(t, (N-1), (M-1));

T_y0 = T_ext * ones(N-1,1);
T_M = 1/3*(4*T(:,end)-T(:,end-1));
T = [T_y0 T T_M];   % Adding boundaries along y
T_N = 1/3*(4*T(end,:)-T(end-1,:));
T_x0 = 1/3*(4*T(1,:)-T(2,:));
T = [T_x0;T;T_N];  %Adding boundaries along x


figure
mesh(y,x,T)
xlabel("y")
ylabel("x")
zlabel("Temperature in metal block, T")
title("Numerical solution when N=120")

   