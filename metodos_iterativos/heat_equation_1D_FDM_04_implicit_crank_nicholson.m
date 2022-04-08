%% Heat equation in 1D

% The PDE for the 1D heat equation is

% u_t = \kappa u_xx, where 0 <= x <= L and t > 0

% where \kappa is the 'heat diffusion' related to the material/concentration

% with boundary conditions
% u(0, t) = a(t)
% u(L, t) = b(t)

% and initial conditions
% u(x, 0) = g(x)

%     x--->
% |---------------------------------------------------------|
% 0                                                         L

%% Discretised 1D heat equation
% The discretised (space and time) version with Finite Differences is

% [U^(i+1) - U^(i)]/dt = \kappa/(2*dx^2) ([U_(j-1)^(i+1) - 2 U_j^(i+1) + U_(j+1)^(i+1)] + [U_(j-1)^(i) - 2 U_j^(i) + U_(j+1)^(i)])

% observe this is an implicit scheme where we are taking 0.5 from U^(i+1)
% values and 0.5 from U^(i) values

% in the domain 0=<x=<L and t >= 0

% Boundary conditions
% U(0,t)=a(t)
% U(L,t)=b(t)

% Initial conditions U(x,0)=g(x)

% We are computing the numerical u(t,x) which is the solution function to
% the heat equation

%% Clear all variables
clear; clc; clf;

%% Problem parameters

% Heat constant
kappa = 0.1;

% Boundary conditions
%a=@(t) 2*t;
%b=@(t) 2*t+1;
a=@(t) 100;
b=@(t) 40;

% Initial condition
%g=@(x) x.^2;
g=@(x) 0;
%g=@(x) sin(x)*cos(x);

% Size of the domain
%L = pi;
L = 1;
% Number of (interior) points (equations) to discretize the domain
Nx = 10;
dx = L/(Nx+1);

% Time integration parameters
dt = 0.1;
initial_time = 0.0;
final_time = 1.0;

% This is no longer needed since we are unsing an implicit strategy
disp(dt)
disp(dx)
disp(dt/dx^2)
% Check for stability condition
%if (dt / dx^2) >= 0.5
 %   disp("The stability condition is not met")
 %   return
%end

%% Problem setup

% Domain discretisation
x = linspace(0, L, Nx+2);
size(x)
disp(x)

% Get the number of time steps
ntime_steps = (final_time - initial_time) / dt;
disp(ntime_steps)

% Create the u solution vector
u=zeros(Nx+2, floor(ntime_steps+1));
% Set initial conditions
for i=1:Nx+2
    u(i, 1) = g(x(i));
end

% Set boundary conditions
u(1, 1) = a(initial_time);
u(Nx+2, 1) = b(initial_time);

%% The factor r
r = kappa * dt/dx^2;

%% Create the matrix
A = eye(Nx);
A = A * (1 + r);
A = A + (-r * 0.5 * diag(ones(Nx-1, 1), -1));
A = A + (-r * 0.5 * diag(ones(Nx-1, 1), 1));

%% Plot the initial conditions
plot(x,u(:,1),'linewidth',2);
hold on

%% Iterate over time
for i = 1:ntime_steps

    % Set boundary conditions for next time step
    u(1,i+1) = u(1,i);
    u(Nx+2,i+1) = u(Nx+2,i);

    % Create an additional vector to keep track of boundary conditions (we
    % should update this vector every iteration if there are time
    % dependencies, in this case there is no need to do so, however we
    % leave it here for future cases)
    gamma=zeros(Nx,1);
    gamma(1) = r * 0.5 * u(1,i) + r * 0.5 * u(1,i+1);
    gamma(Nx) = r * 0.5 * u(Nx+2,i) + r * 0.5 * u(Nx+2,i+1);

    % The RHS
    rhs = zeros(Nx,1);
    rhs(1) = (1 - r) * u(2, i) + r * 0.5 * u(3, i);
    for j = 2:Nx-1
        rhs(j) = r * 0.5 * u(j, i) + (1 - r) * u(j+1, i) + r * 0.5 * u(j+2, i);
    end
    rhs(Nx) = r * 0.5 * u(Nx, i) + (1 - r) * u(Nx+1, i);

    % Add boundary conditions
    rhs = rhs + gamma;

    % Solve the system of equations
    sol = linsolve(A, rhs);
    %sol = A\rhs;

    % This is for ploting ------------------------------
    % Copy the solution into the vector u
    u(2:Nx+1, i+1) = sol;

    % Add boundary conditions
    u(1,i+1) = u(1,i);
    u(Nx+2,i+1) = u(Nx+2,i);
    
    plot(x,u(:,i+1),'linewidth',2);
    %plot(x(2:Nx)',u,'linewidth',2);
    hold on

end

%% Plot results
title('Heat Equation 1D (views for different times)')
xlabel('x')
ylabel('u')
%times = initial_time:dt:final_time;
%legend(num2str(times'),'Location','northwest');

