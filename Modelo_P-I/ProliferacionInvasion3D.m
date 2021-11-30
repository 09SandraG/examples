function [u,x] = ProliferacionInvasion3D(m,n,p,t,f,v,a,b,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretizacion mediante el método de diferencias finitas %%
% Para el modelo de proliferacion invasion.                %%
% Sandra I. García Mendoza                                 %%
% Noviembre 2021                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Se hace un metodo iterativo para la aproximacion.
%
% Argumentos de entrada
%   m     Entero       Numero de nodos en la malla en el espacio en x.
%   n     Entero       Numero de nodos en la malla en el espacio en y.
%   p     Entero       Numero de nodos en la malla en el espacio en z.
%   t     Entero       Numero de pasos en el tiempo
%   v     Real         Coeficiente de difusion
%   a     Real         Velocidad de propagacion en el eje x.
%   b     Real         Velocidad de propagacion en el eje y.
%   c     Real         Velocidad de propagacion en el eje z.
%   f     funcion      funcion para condicion inicial y de frontera
%
% Argumentos de salida
%   u_approx     matriz       Matriz con la solucion aproximada
%   x            matriz       Vector de coordenadas x de la malla
%   y            matriz       Vector de coordenadas y de la malla
%   z            matriz       Vector de coordenadas z de la malla
% 
% Ejemplo de uso:
%       [u,x] = ProliferacionInvasion(21,21,21,1000,@fAdvDif,0.2,0.3,0.3,0.3);
%
%
% Se inicializan las variables
close all                            % Cerramos todas las ventanas de figuras.
T = linspace(0,1,t);                % Se crea la malla temporal.
x = linspace(0,1,m);              % Se crea la malla espacial.
y = linspace(0,1,n);              % Se crea la malla espacial.
z = linspace(0,1,p);              % Se crea la malla espacial.
[x,y,z] = meshgrid(x,y,z);
dt = T(2) - T(1);                    % Se calcula dt.
dx = x(2,1,1) - x(1,1,1);                    % Se calcula dx.
dy = y(1,2,1) - y(1,1,1);                    % Se calcula dx.
dz = z(1,1,2) - z(1,1,1);                    % Se calcula dx.
u = zeros(m,n,p,t);               % Se inicializa u con ceros.
ur = zeros(m,n,p,t);               % Se inicializa ur con ceros.
%c = abs(a*dt/dx + b*dt/dy + c*dt/dz);                     % Calculamos c para ver la estabilidad.
%if c > 1/2
%    fprintf('El codigo puede no funcionar c = %f\n',c);
%end

%Agregamos condiciones inicial
for i = 1:m
    for j = 1:n
        for k = 1:p
            u(i,j,k,1) = f(x(i,j,k),y(i,j,k),z(i,j,k),T(1),v,a,b,c);  % Se agrega la condicion inicial.
        end
    end
end

% Agregamos la condicion de frontera.
for k = 1:t
    for j = 1:n
        u(i,1,k) = f(x(i,1),y(i,1),T(k),v,a,b);  % Se agrega la condicion de frontera en y=0.
        u(i,n,k) = f(x(i,n),y(i,n),T(k),v,a,b);  % Se agrega la condicion de frontera en y=1.
    end
    for j = 2:n
        u(i,1,k) = f(x(i,1),y(i,1),T(k),v,a,b);  % Se agrega la condicion de frontera en x=0.
        u(i,n,k) = f(x(i,n),y(i,n),T(k),v,a,b);  % Se agrega la condicion de frontera en x=1.
    end
end

% Hacemos el metodo para Diferencias de segundo orden
for k = 2:t
    for i = 2:m-1
        for j = 2:n-1
            u(i,j,k) = u(i,j,k-1) + (v*dt/(dx^2))*(u(i+1,j,k-1) - 2*u(i,j,k-1) + u(i-1,j,k-1)) + (v*dt/(dy^2))*(u(i,j+1,k-1) - 2*u(i,j,k-1) + u(i,j-1,k-1)) - (a*dt/(2*dx))*(u(i+1,j,k-1) - u(i-1,j,k-1) - (b*dt/(2*dy))*(u(i,j+1,k-1)) - u(i,j-1,k-1)); % Se hace el llenado de u para todos los puntos.
        end
    end
end

% Calculamos la solucion exacta
for k = 1:t
    for i = 1:m
        for j = 1:n
            ur(i,j,k) = f(x(i,j),y(i,j),T(k),v,a,b);
        end
    end
end

% Graficamos la solucion 
scrsz = get(groot, 'ScreenSize');  % Se obtienen los limites de la
figure('OuterPosition',[1 1 scrsz(4)]); % Se crea una figura del tamaño
p = ceil(t/100);                     % Se calcula p para no graficar todos los 
for i = 1:p:t
    subplot(1,2,1)
    surf(x,y,u(:,:,i));           % Se grafican la solucion aproximada.
    title('Aproximacion')
    axis([0 1 0 1 -1 1])          % Se ajustan los limites de la grafica.
    subplot(1,2,2)
    surf(x,y,ur(:,:,i));          % Se grafican la solucion exacta 
    title('solucion exacta')
    axis([0 1 0 1 -1 1])          % Se ajustan los limites de la grafica.
    pause(0.01)
end
end