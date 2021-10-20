function [phi_approx,phi_exacta,x,y,z] = Poisson3D2(m1,m2,m3,phi,f)
% funcion que calcula una aproximacion a la solucon de la ecuacion de
% Poisson en 2D.
%
% phi_xx + phi_yy + phi_zz = -f
%
% Se hace un método iterativo para la aproximacion
%
% Escuela Nacional de Optimizacion y Analisis Numerico 2021
%
% Argumentos de entrada
%   m     Entero       Numero de nodos en la malla
%   phi   Funcion      Funcion para condiciones de frontera.
%   f     Funcion      Lado derecho de la ecuacion de Poisson.
%
% Argumentos de salida
%   phi_approx  Matriz  Matriz con la solucion aproximada.
%   x           Matriz  Matriz de coordenadas x de la malla.
%   y           Matriz  Matriz de coordenadas y de la malla.
%
% Ejemplo de uso:
%          [phi_approx,phi_exacta,x,y,z] = Poisson3D2(11,11,11,@phi,@f);
% Inicializamos las variables
close all
m1=11; m2=11; m3=11;
x = linspace(0,1,m1)';                  % Se crea la discretizacion en x.
y = linspace(0,1,m2)';                  % Se crea la discretizacion en y.
z = linspace(0,1,m3)';                  % Se crea la discretizacion en z.
h = x(2) - x(1);                      % Se calcula h.
[x,y,z] = meshgrid(x,y,z);                % Creamos la malla completa. 
phi_approx = zeros(m1,m2,m3);              % Se inicializa la phi_approx.
err = 1;                              % Se incializa err con 1.
tol = sqrt(eps);                      % Se impone una tolerancia para err.
%tol = 0.0006;

% Agregamos condiciones de frontera
for i = 1:m1
    phi_approx(i,1,1) = phi(x(i,1,1),y(i,1,1),z(i,1,1)); % Se agrega la condicion de forntera inferior.
    phi_approx(i,m2,1) = phi(x(i,m2,1),y(i,m2,1),z(i,m2,1)); % Se agrega la condicion de forntera superior.
    phi_approx(i,1,m3) = phi(x(i,1,m3),y(i,1,m3),z(i,1,m3)); % Se agrega la condicion de forntera inferior.
    phi_approx(i,m2,m3) = phi(x(i,m2,m3),y(i,m2,m3),z(i,m2,m3)); % Se agrega la condicion de forntera superior.
end
for i = 1:m2
    phi_approx(1,i,1) = phi(x(1,i,1),y(1,i,1),z(1,i,1)); % Se agrega la condicion de forntera izquierda.
    phi_approx(m1,i,1) = phi(x(m1,i,1),y(m1,i,1),z(m1,i,1)); % Se agrega la condicion de forntera derecha.
    phi_approx(1,i,m3) = phi(x(1,i,m3),y(1,i,m3),z(1,i,m3)); % Se agrega la condicion de forntera izquierda.
    phi_approx(m1,i,m3) = phi(x(m1,i,m3),y(m1,i,m3),z(m1,i,m3)); % Se agrega la condicion de forntera derecha.
end
for i = 1:m3
    phi_approx(1,1,i) = phi(x(1,1,i),y(1,1,i),z(1,1,i)); % Se agrega la condicion de forntera izquierda.
    phi_approx(m1,1,i) = phi(x(m1,1,i),y(m1,1,i),z(m1,1,i)); % Se agrega la condicion de forntera derecha.
    phi_approx(1,m2,i) = phi(x(1,m2,i),y(1,m2,i),z(1,m2,i)); % Se agrega la condicion de forntera izquierda.
    phi_approx(m1,m2,i) = phi(x(m1,m2,i),y(m1,m2,i),z(m1,m2,i)); % Se agrega la condicion de forntera derecha.
end


% Llenamos la matriz usando Diferencias Finitas
while err >= tol
    err = 0;
    for i = 2:m1-1
        for j = 2:m2-1                   % Se utilizan Diferencias centradas
            for k = 2:m3-1
                t = (1/6)*(phi_approx(i-1,j,k) + phi_approx(i+1,j,k) + phi_approx(i,j-1,k) + phi_approx(i,j+1,k) + phi_approx(i,j,k-1) + phi_approx(i,j,k+1) - h^2*f(x(i,j,k),y(i,j,k),z(i,j,k)));
                err = max(err,abs(t-phi_approx(i,j,k)));% Se calcula un nuevo err.
                phi_approx(i,j,k) = t;           % Se asigna el valor de phi_approx.
            end
        end
    end
end

% Calculamos la solucion exacta
phi_exacta = phi(x,y,z);

% GRaficamos la solucion
% scrsz = get(groot,'ScreenSiza');          % Se obtienen los limites de la
% figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]); % Se crea una figura del
% tamaño de 
%figure                                 
%subplot(1,2,1)                        % Se divide la grafica en 2.
%surf(x,y,phi_approx(:,:,1));                 % Se grafica la solucion aproximada.
%title('Aproximacion')                 % Se agrega el titulo de la grafica.
%subplot(1,2,2)                        % Se usa la otra parte de la grafica.
%surf(x,y,phi_exacta);                 % Se grafica la solucion exacta.
%title('Solucion exacta')              % Se agrega titulo a la grafica.
end