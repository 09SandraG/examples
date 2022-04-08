%%Modelo de proliferación-invasión 2D

% c_t = Dif*(c_xx + c_yy) + g
% donde 0 <= x <= M, 0 <= y <= N, t>0,
% Dif es el coeficiente de difusión,
% g denota la función de crecimiento y tratamiento (en caso de tener
% alguno)

% Condiciones de frontera 
% c(0,y,t) = a1(t);  %0
% c(M,y,t) = b1(t);  %0
% c(x,0,t) = a2(t);  %0
% c(x,N,t) = b2(t);  %0

% Condiciones de inicio
% c(x,y,0) = h(x,y);  %4000

%% Discretización de la ecuación en 2D

% c^(k+1)_(i,j) = c^(k)_(i,j) + Dif*dt*(c^(k)_(i+1,j) - 2*c^(k)_(i,j) + c^(k)_(i-1,j))/dx^2 +
% Dif*dt*(c^(k)_(i,j+1) - 2*c^(k)_(i,j) + c^(k)_(i,j-1))/dy^2 + g^(k)_(i,j)

%% Clear all variables
clear; clc; clf;

%% Inicialización de los parámetros

dif = 0.003; % Coeficiente de difusión
r = 0.012;
g = @(c) 0; % r*c;

%Condiones de frontera
a1 = @(t) 0;
b1 = @(t) 0;
a2 = @(t) 0;
b2 = @(t) 0;

%Condiciones de inicio
h = @(x,y) 4000;

%Parámetros de discretización sobre el dominio del espacio
M = 200;  %Tamaño del dominio sobre x
N = 200;  %Tamaño del dominio sobre y
mx = 30;  %Número de puntos para la discretización sobre x
ny = 30;  %Número de puntos para la discretización sobre y
dx = M/(mx+1);   %Distancia entre los puntos interiores sobre x
dy = N/(ny+1);   %Distancia entre los puntos interiores sobre y

%Parámetros de discretización sobre el dominio del tiempo
T = 10;   %Tiempo total transcurrido (Número de días)
dt = 0.01; %Tamaño de paso de tiempo;

%Condición de estabilidad
if dt/dx >= 1/sqrt(2)     %para dx=dy     %dt/dx >= sqrt(3/8); 
    fprintf('No se cumple la condión de estabilidad');
    return
end

%% Desarrollo del problema

% Discretización del dominio espacial
x = linspace(0,M,mx+1);              % Se crea la malla espacial en x.
y = linspace(0,N,ny+1);              % Se crea la malla espacial en y.
[x,y] = meshgrid(x,y);             % Se cre la malla espacial

% Condiciones de inicio
c = zeros(mx+1,ny+1);
for i=5:25 %2:mx
    for j=5:25 %2:ny
        c(i,j) = h(x(i,j),y(i,j));
    end 
end

% Condiciones de frontera
for j=1:ny+1
    c(1,j) = a1(0);
    c(mx+1,j) = b1(0);
end
for i=1:mx+1
    c(i,1) = a2(0);
    c(i,ny+1) = b2(0);
end

%% Grafica del problema con CI
posPlot = 1;
subplot(1,6,posPlot)
surf(x,y,c) 
xlabel('x')
ylabel('y')
zlabel('z')
title('c en el tiempo inicial')

%% Iteraciones sobre el tiempo
c_approx = zeros(mx+1,ny+1);
for t=0:dt:T
    
    for i=2:mx
        for j=2:ny
            c_approx(i,j) = c(i,j) + (dif*dt/dx^2)*(c(i+1,j) - 2*c(i,j) + c(i-1,j)) + (dif*dt/dy^2)*(c(i,j+1) - 2*c(i,j) + c(i,j-1)) + g(c(i,j));
        end 
    end
    
    % Agregamos condiciones de frontera
    for j=1:ny
        c_approx(1,j) = a1(t);
        c_approx(mx+1,j) = b1(t);
    end
    for i=1:mx
        c_approx(i,1) = a2(t);
        c_approx(i,ny+1) = b2(t);
    end
    
    % Graficamos 
    mdt = T/5;
    csp = mod(t,mdt);
    if csp == 0 && t ~= 0 && posPlot < 6
        posPlot = posPlot +1;
        subplot(1,6,posPlot)
        surf(x,y,c_approx,'FaceAlpha',0.5)
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title(['c en el tiempo ', num2str(t)])
    end
    
    % Transferir la solución a c
    c(2:mx,2:ny) = c_approx(2:mx,2:ny);
   
end
    
% Grafica de resultados


