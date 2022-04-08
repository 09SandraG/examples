%% Heat equation  in 1D

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
% The discretised (in space and time) version with Finite Differences is

% [U^(j+1) - U^(j)]/dt = \kappa [U_(i-1) - 2 U_i + U_(i+1)]/dx^2

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
% Number of points (equations) to discretize the domain
Nx = 10;
dx = L/(Nx+1);

% Time integration parameters
dt = 0.003;
initial_time = 0;
final_time = 1;

disp(dx)
disp(dt)
disp(dt/dx^2)
% Check for stability condition
 if (dt / dx^2) >= 0.5
     disp("The stability condition is not met")
     return
 end

%% Problem setup

% Domain discretisation
x = linspace(0, L, Nx+1);

% Set initial conditions
u=zeros(Nx+1,1);
for i=2:Nx
    u(i) = g(x(i));
end

% Set boundary conditions
u(1) = a(initial_time);
u(Nx+1) = b(initial_time);

%% Iterate over time
u_new=zeros(Nx+1,1);
for t=initial_time:dt:final_time

    for i=2:Nx
        u_new(i) = u(i) + dt * (kappa/dx^2 * (u(i-1) - 2*u(i) + u(i+1)));
    end
    
    % Add boundary conditions for plotting
    uplot = zeros(Nx+1,1); % Create the vector for plotting
    uplot(1) = a(t); % BC
    uplot(2:Nx) = u_new(2:Nx); % The solution at t
    uplot(Nx+1) = b(t); % BC
    
    plot(x,uplot,'linewidth',2);
    %plot(x(2:Nx)',u,'linewidth',2);
    hold on

    % Transefer the soltuion to vector u
    u(2:Nx) = u_new(2:Nx);

end

%% Plot results
title('Heat Equation 2D views for different times')
xlabel('x')
ylabel('u')
%times = initial_time:dt:final_time;
%legend(num2str(times'),'Location','northwest');

