function [phi_approx,phi_exacta,x,y,tiempo] = Poisson2D2(m,n,phi,f)
% funcion que calcula una aproximacion a la solucon de la ecuacion de
% Poisson en 2D.
%
% phi_xx + phi_yy = -f
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
%          [phi_approx,phi_exacta,x,y] = Poisson2D2(11,11,@phi,@f);
% Inicializamos las variables
close all
x = linspace(0,1,m)';                  % Se crea la discretizacion en x.
y = linspace(0,1,n)';                  % Se crea la discretizacion en y.
h = x(2) - x(1);                      % Se calcula h.
[x,y] = meshgrid(x,y);                % Creamos la malla completa. 
phi_approx = zeros(m,n);              % Se inicializa la phi_approx.
err = 1;                              % Se incializa err con 1.
tol = sqrt(eps);                      % Se impone una tolerancia para err.
%tol = 0.0006;

% Agregamos condiciones de frontera
for i = 1:m
    phi_approx(i,1) = phi(x(i,1),y(i,1)); % Se agrega la condicion de forntera inferior.
    phi_approx(i,n) = phi(x(i,n),y(i,n)); % Se agrega la condicion de forntera superior.
end
for i = 1:n
    phi_approx(1,i) = phi(x(1,i),y(1,i)); % Se agrega la condicion de forntera izquierda.
    phi_approx(m,i) = phi(x(m,i),y(m,i)); % Se agrega la condicion de forntera derecha.
end
cont11 = 1;
tStart = cputime;
% Llenamos la matriz usando Diferencias Finitas
while (err >= tol) && (cont11 < 1500)
    err = 0;
    for i = 2:m-1
        for j = 2:n-1                   % Se utilizan Diferencias centradas
            t = (1/4)*(phi_approx(i-1,j) + phi_approx(i+1,j) + phi_approx(i,j-1) + phi_approx(i,j+1) - h^2*f(x(i,j),y(i,j)));
            err = max(err,abs(t-phi_approx(i,j)));% Se calcula un nuevo err.
            phi_approx(i,j) = t;           % Se asigna el valor de phi_approx.
        end
    end
    cont11=cont11+1;
end
tiempo = cputime - tStart;


% Calculamos la solucion exacta
phi_exacta = phi(x,y);
fprintf('Despues de %3.0f iteraciones el error de la aproximación es: %3.6e\n',cont11,err)
% Graficamos la solucion
% scrsz = get(groot,'ScreenSiza');          % Se obtienen los limites de la
% figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]); % Se crea una figura del
% tamaño de 
%figure                                 
%subplot(1,2,1)                        % Se divide la grafica en 2.
%surf(x,y,phi_approx);                 % Se grafica la solucion aproximada.
%title('Aproximacion')                 % Se agrega el titulo de la grafica.
%subplot(1,2,2)                        % Se usa la otra parte de la grafica.
%surf(x,y,phi_exacta);                 % Se grafica la solucion exacta.
%title('Solucion exacta')              % Se agrega titulo a la grafica.
end