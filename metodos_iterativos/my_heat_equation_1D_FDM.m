%% Heat equation  in 1D

% The PDE for the 1D heat equation is

% u_t = \kappa u_xx, where 0 <= x <= L and t > 0
% \kappa is a 'heat' constant

% with boundary conditions
% u(0, t) = a(t)
% u(L, t) = b(t)

% and initial conditions
% u(x, 0) = g(x)

%     x--->
% |---------------------------------------------------------|
% 0                                                         L

%% Discretised 1D heat equation
% The discretised (space) version with Finite Differences is

% dU_i/dt = \kappa [U_(i-1) - 2 U_i + U_(i+1)]/dx^2

% The PDE for 1D heat equation is Ut = \kappa Uxx
% in the domain 0=<x=<L and 0=<t

% Boundary conditions
% U(0,t)=a(t)
% U(L,t)=b(t)

% Initial conditions U(x,0)=g(x)

% u(t,x) is the solution

%% Clear all variables
clear; clc; clf;

%% Problem configuration

% Heat constant
kappa = 0.01;
%a=@(t) 2*t;
%b=@(t) 2*t+1;
a=@(t) 0;
b=@(t) 0;
%g=@(x) x.^2;
%g=@(x) sin(x)*cos(x);
g=@(x) 2;

% Size of the domain
L = pi;
% Number of points (equations) to discretize the domain
Nx = 100;
dx = L/(Nx-1);

% Domain discretisation
x = linspace(0, 1, Nx);

% Time integration parameters
dt = 0.1;
initial_time = 0;
final_time = 1;

% Set initial conditions
u=zeros(Nx,1);
for i=2:Nx-1
    u(i) = g(x(i));
end

% Set boundary conditions
u(1) = a(initial_time);
u(Nx) = b(initial_time);

%% Create the system of equations
% Matrix A
A = eye(Nx);
A = A * (-2 / dx^2);
A = A + 1 / dx^2 * diag(ones(Nx-1, 1), -1);
A = A + 1 / dx^2 * diag(ones(Nx-1, 1), 1);
A = kappa * A;

for t=initial_time:dt:final_time
    % Vector rhs (from boundary and initial conditions)
    rhs = -(kappa * u);

    % Solve the system
    du = A\rhs;

    % Time integration (Euler)
    u = u + dt*du;

    % Enforce boundary conditions
    %u(1) = a(t);
    %u(Nx) = b(t);
    
    plot(x,u,'linewidth',2);
    hold on

end

title('Heat Equation 2D views for different times')
xlabel('x')
ylabel('u')
times = initial_time:dt:final_time;
legend(num2str(times'),'Location','northwest');