function [tiempo_sum,iter,error] = ProliferacionInvasion2D_MetodoImplicito(m,n)
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

% (1/dt)*(c^(k+1)_(i,j) - c^(k)_(i,j)) = Dif*(1/(2*dx^2))*(c^(k+1)_(i+1,j) - 2*c^(k+1)_(i,j) + c^(k+1)_(i-1,j)) +
% Dif*(1/(2*dx^2))*(c^(k)_(i+1,j) - 2*c^(k)_(i,j) + c^(k)_(i-1,j)) + 
% Dif*(1/(2*dy^2))*(c^(k+1)_(i,j+1) - 2*c^(k+1)_(i,j) + c^(k+1)_(i,j-1)) + 
% Dif*(1/(2*dy^2))*(c^(k)_(i,j+1) - 2*c^(k)_(i,j) + c^(k)_(i,j-1)) + g^(k)_(i,j)

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
mx = m; %30;  %Número de puntos para la discretización sobre x
ny = n; %30;  %Número de puntos para la discretización sobre y
dx = M/(mx+1);   %Distancia entre los puntos interiores sobre x
dy = N/(ny+1);   %Distancia entre los puntos interiores sobre y

%Parámetros de discretización sobre el dominio del tiempo
T = 10;   %Tiempo total transcurrido (Número de días)
dt = 0.1; %Tamaño de paso de tiempo;

%Parámetros para la discretización
% lambda se considerará con el valor de 1/2
rx = dif*dt/dx^2;
ry = dif*dt/dy^2;

%% Desarrollo del problema

% Discretización del dominio espacial
x = linspace(0,M,mx+1);              % Se crea la malla espacial en x.
y = linspace(0,N,ny+1);              % Se crea la malla espacial en y.
[x,y] = meshgrid(x,y);             % Se cre la malla espacial

% Condiciones de inicio  
c = zeros(mx-1,ny-1);
for i=mx/5:(4*mx)/5   %2:mx-1
    for j=(ny)/5:(4*ny)/5  %2:ny-1
        c(i,j) = h(x(i,j),y(i,j));
    end
end

% Condiciones de frontera (Se guardan en b)
% Para guardar las condiciones de frontera sobre x
bcx = zeros((mx-1)*(ny-1),1); 
for j=1:ny-1
    bcx(j) = rx*a1(0);
    bcx((mx-1)*(ny-1)-(ny-1)+j) = rx*b1(0);
end
%disp(bcx)

% Para guardar las condiciones de frontera sobre y
bcy = zeros((mx-1)*(ny-1),1);
for i=1:ny-1:(mx-1)*(ny-1)
    bcy(i) = ry*a2(0);
    bcy(i+(ny-2)) = ry*b2(0);
end
%disp(bcy)

%Para guardar la evaluación de g  
vg = zeros((mx-1)*(ny-1),1);
for i=0:(mx-2)
    for j=1:(ny-1)
        vg(i*(ny-1)+j) = dt*g(c(i+1,j));
    end
end

% Suma de los vectores
bc = bcx + bcy + vg;

%% Grafica del problema con CI 

% En c_plot agregamos todos los valores que se van a graficar
c_plot = zeros(mx+1,ny+1);
c_plot(1,:) = rx*a1(0);
c_plot(:,1) = ry*a2(0);
c_plot(2:mx,2:ny) = c;
c_plot(mx,:) = rx*b1(0);
c_plot(:,ny) = ry*b2(0);

% Grafica en el tiempo 0;
posPlot = 1;
%%%%%%%%%%subplot(1,6,posPlot)
%figure
%heatmap(c_plot)
%%%%%%%%%%%%%surf(x,y,c_plot) 
%xlabel('x')
%ylabel('y')
%%%%%%%%%%%%zlabel('z')
%title('c en el tiempo inicial')

%% Creación de matrices
% MATRIZ PARA EL LADO DERECHO (Donde se hace evaluación)
% Para x
Ax_d = eye((mx-1)*(ny-1));
Ax_d = Ax_d * (-rx);
Ax_d = Ax_d + ((rx/2)*diag(ones((mx-1)*(ny-1)-ny, 1), -ny));
Ax_d = Ax_d + ((rx/2)*diag(ones((mx-1)*(ny-1)-ny, 1), ny));
%spy(Ax)

%Para y
Ay_d = zeros((mx-1)*(ny-1));
Ayp_d = eye(ny-1);
Ayp_d = Ayp_d * (-ry);
Ayp_d = Ayp_d + ((ry/2)*diag(ones(ny-2, 1), -1));
Ayp_d = Ayp_d + ((ry/2)*diag(ones(ny-2, 1), 1));
j=1;
for i=1:ny-1:(mx-1)*(ny-1)
    Ay_d(i:j*(ny-1),i:j*(ny-1)) = Ayp_d;
    j=j+1;
end
%spy(Ay)

A_der = eye((mx-1)*(ny-1));
A_der = A_der + Ax_d + Ay_d;

% MATRIZ PARA EL LADO IZQUIERDO (Donde se va a resolver)
% Para x
Ax_i = eye((mx-1)*(ny-1));
Ax_i = Ax_i * (rx);
Ax_i = Ax_i + ((-rx/2)*diag(ones((mx-1)*(ny-1)-ny, 1), -ny));
Ax_i = Ax_i + ((-rx/2)*diag(ones((mx-1)*(ny-1)-ny, 1), ny));
%spy(Ax)

%Para y
Ay_i = zeros((mx-1)*(ny-1));
Ayp_i = eye(ny-1);
Ayp_i = Ayp_i * (ry);
Ayp_i = Ayp_i + ((-ry/2)*diag(ones(ny-2, 1), -1));
Ayp_i = Ayp_i + ((-ry/2)*diag(ones(ny-2, 1), 1));
j=1;
for i=1:ny-1:(mx-1)*(ny-1)
    Ay_i(i:j*(ny-1),i:j*(ny-1)) = Ayp_i;
    j=j+1;
end
%spy(Ay)

A_izq = eye((mx-1)*(ny-1));
A_izq = A_izq + Ax_i + Ay_i;


%% Iteraciones sobre el tiempo
c_approx = zeros(mx+1,ny+1);
vector_b = zeros((mx-1)*(ny-1),1);
tiempo = zeros(T*10,1);
iteraTotal = zeros(T*10,1);
err_Total = zeros(T*10,1);
t_aux = 0;
for t=0:dt:T
    % Variable auxiliar para toma de tiempos
    t_aux = t_aux + 1;
    
    % Cambiamos de matriz c a vector vc
    vc = reshape(c,(mx-1)*(ny-1),1); 
    
    % Evaluación del lado derecho del sistema algebraico
    vector_b = A_der*vc + bc;
    
    % Resolver A*x=b
    tStart = cputime;
    %vc=linsolve(A_izq,vector_b);    % Método definido por Matlab
    %Método LU
    [L,U,P] = lu(A_izq);           % Factorizamos la matriz en una triangular inferior y una superior
    vM = L\(P*vector_b);                  % Resolvemos la matriz inferior
    vc = U\vM;
    %Método iterativo. Gauss-Seidel
    %addpath('C:\Users\Dell\Documents\Curso_Computo_Cientifico\Poisson\examples\metodos_iterativos');
    %[vc,num_iter,err] = GaussSeidel(A_izq,vector_b,vc,(mx-1)*(ny-1),1e-14,500);  %% Método implementado por mí
    %vc = Gauss_Siedel(A_izq,vector_b,vc,0.000001);   %% Implementación descargada de internet(https://www.mathworks.com/matlabcentral/fileexchange/73488-gauss-seidel-iterative-method)
    %[vc,fl,err,num_iter,rv] = pcg(A_izq,vector_b,1e-14,500);
    tiempo(t_aux) = cputime - tStart;
    %iteraTotal(t_aux) = num_iter;
    %err_Total(t_aux) = err;
    
   
    % Agregamos condiciones de frontera
    ctemp = reshape(vc,mx-1,ny-1); 
    for j=1:ny+1
        c_approx(1,j) = rx*a1(t);
        c_approx(mx+1,j) = rx*b1(t);
    end
    for i=1:mx+1
        c_approx(i,1) = ry*a2(t);
        c_approx(i,ny+1) = ry*b2(t);
    end
    c_approx(2:mx,2:ny) = ctemp;
    
    % Graficamos 
    %mdt = T/5;
    %csp = mod(t,mdt);
    %if csp == 0 && t ~= 0 && posPlot < 6
        %posPlot = posPlot +1;
        %%%%subplot(1,6,posPlot)
        %figure
        %heatmap(c_approx)
        %%%%surf(x,y,c_approx,'FaceAlpha',0.5)
        %xlabel('x')
        %ylabel('y')
        %%%%zlabel('z')
        %title(['c en el tiempo ', num2str(t)])
    %end
    
    c = c_approx(2:mx,2:ny);
   
end

%tiempo_prom = mean(tiempo);    % Para hacer comparaciones de tiempos
%promedios don diferentes métodos numéricos
tiempo_sum = sum(tiempo);       % Para hacer comparaciones de tiempos 
%totales entre método explicito e implicito
iter = mode(iteraTotal);
error = mean(err_Total);
end


