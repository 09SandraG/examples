function [tiempo_sum,tiempo_Prom,tiempo_moda,rr] = ProliferacionInvasion2D_MetodoExplicito_ConMatrices(mx,ny)
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
%clear; clc; clf;

%% Inicialización de los parámetros

dif = 0.003; % Coeficiente de difusión
r = 0.012;
g = @(c) r*c;

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
%mx = 30;  %Número de puntos para la discretización sobre x
%ny = 30;  %Número de puntos para la discretización sobre y
dx = M/(mx+1);   %Distancia entre los puntos interiores sobre x
dy = N/(ny+1);   %Distancia entre los puntos interiores sobre y

%Parámetros de discretización sobre el dominio del tiempo
T = 10;   %Tiempo total transcurrido (Número de días)
dt = 0.001; %Tamaño de paso de tiempo;

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
c = zeros(mx-1,ny-1);
for i=mx/5:(4*mx)/5   %2:mx-1
    for j=(ny)/5:(4*ny)/5   %j=2:ny-1
        c(i,j) = h(x(i,j),y(i,j));
    end
end


% Condiciones de frontera (Se guardan en b)
% Para guardar las condiciones de frontera sobre x
bcx = zeros((mx-1)*(ny-1),1); 
for j=1:ny-1
    bcx(j) = a1(0);
    bcx((mx-1)*(ny-1)-(ny-1)+j) = b1(0);
end
%disp(bcx)

% Para guardar las condiciones de frontera sobre y
bcy = zeros((mx-1)*(ny-1),1);
for i=1:ny-1:(mx-1)*(ny-1)
    bcy(i) = a2(0);
    bcy(i+(ny-2)) = b2(0);
end
%disp(bcy)

%Para guardar la evaluación de g  %%%%%%%REVISAR%%%%%%%%%%%
vg = zeros((mx-1)*(ny-1),1);
for i=0:(mx-2)
    for j=1:(ny-1)
        vg(i*(ny-1)+j) = g(c(i+1,j));
    end
end

% Suma de los vectores
bc = bcx + bcy + vg;

%% Grafica del problema con CI 

% En c_plot agregamos todos los valores que se van a graficar
c_plot = zeros(mx+1,ny+1);
c_plot(1,:) = a1(0);
c_plot(:,1) = a2(0);
c_plot(2:mx,2:ny) = c;
c_plot(mx,:) = b1(0);
c_plot(:,ny) = b2(0);

% Grafica en el tiempo 0;
posPlot = 1;
%subplot(1,6,posPlot)
%surf(x,y,c_plot) 
%xlabel('x')
%ylabel('y')
%zlabel('z')
%title('c en el tiempo inicial')

%% Creación de matrices
% Para x

Ax = eye((mx-1)*(ny-1));
Ax = Ax * (-2);
Ax = Ax + diag(ones((mx-1)*(ny-1)-ny, 1), -ny);
Ax = Ax + diag(ones((mx-1)*(ny-1)-ny, 1), ny);
%spy(Ax)

%Para y
Ay = zeros((mx-1)*(ny-1));
Ayp = eye(ny-1);
Ayp = Ayp * (-2);
Ayp = Ayp + diag(ones(ny-2, 1), -1);
Ayp = Ayp + diag(ones(ny-2, 1), 1);
j=1;
for i=1:ny-1:(mx-1)*(ny-1)
    Ay(i:j*(ny-1),i:j*(ny-1)) = Ayp;
    j=j+1;
end
%spy(Ay)

%A = Ax+Ay;
vc = reshape(c,(mx-1)*(ny-1),1); % Cambiamos de matriz a vector
%% Iteraciones sobre el tiempo
tiempo = zeros(T*1000,1);
ti_aux=0;
c_approx = zeros(mx+1,ny+1);
dudt = zeros((mx-1)*(ny-1),1);
for t=0:dt:T
    
    
    ti_aux = ti_aux + 1;
    tStart = cputime;            % Toma de tiempos
    % Evaluación de las derivadas en el tiempo t
    %dudt = (dif/dx^2)*Ax*vc + (dif/dy^2)*Ay*vc + bc;
    
    % Integración en el tiempo
    %vc = vc + dt*dudt;
    vc = vc + dt*((dif/dx^2)*Ax*vc + (dif/dy^2)*Ay*vc + bc);
    tiempo(ti_aux) = cputime - tStart; % guardamos los tiempos
    %fprintf('%6.3f %6.2f %12.12f\r\n', t, ti_aux, tiempo(ti_aux));
    
    
    % Agregamos condiciones de frontera
    ctemp = reshape(vc,mx-1,ny-1); 
    for j=1:ny+1
        c_approx(1,j) = a1(t);
        c_approx(mx+1,j) = b1(t);
    end
    for i=1:mx+1
        c_approx(i,1) = a2(t);
        c_approx(i,ny+1) = b2(t);
    end
    c_approx(2:mx,2:ny) = ctemp;
    
    % Graficamos 
    %mdt = T/5;
    %csp = mod(t,mdt);
    %if csp == 0 && t ~= 0 && posPlot < 6
    %    posPlot = posPlot +1;
    %    subplot(1,6,posPlot)
    %    surf(x,y,c_approx,'FaceAlpha',0.5)
    %    xlabel('x')
    %    ylabel('y')
    %    zlabel('z')
    %    title(['c en el tiempo ', num2str(t)])
    %end
    
   
end
   
%fprintf('%12.8f\n', tiempo');
tiempo_sum = sum(tiempo);       % Para hacer comparaciones de tiempos 
tiempo_Prom = mean(tiempo);
tiempo_moda = mode(tiempo);
rr = unique(tiempo);
% Grafica de resultados


