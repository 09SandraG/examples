function [phi_approx, phi_exacta,x,y,tiempo] = Poisson2D(m,n,phi,f)
% Función que calcula una aproximación a la solución de la ecuación de
% Possion en 2D.
%
%phi_xx + phi_yy = -f
%
%Se ensamblan las matrices de Diferencias para resolver el sistema lineal.
%
%Escuela Nacional de Optimización y Análisis Numérico 2021
%
%Argumentos de entrada
%   m     Entero     Número de nodos en la malla en x.
%   phi   función    Función para condiciones de frontera.
%   f     función    Lado derecho de la ecuación de Poisson.
%
%Argumentos de salida 
%   phi_approx   Matriz    Matriz de la solución aproximada.
%   x            Matriz    Matriz de coordenadas x de la malla.
%   y            Matriz    Matriz de coordenadas y de la malla.
%
%Ejemplo de uso:
% [phi_aprox,phi_exacta,x,y] = Poisson2D(11,11,@phi,@f);
%
% Inicializamos las variables
close all                          % Cierra ventanas de figuras abiertas
x = linspace(0,1,m);              % Se crea la discretización en x
y = linspace(0,1,n);              % Se crea la discretización en y
h = x(2) - x(1);                   % Se calcula h
[x,y] = meshgrid(x,y);             % Creamos la malla completa.
A = zeros((m-2)*(n-2),(m-2)*(n-2));% Se inicializa A con ceros.
rhs = zeros((m-2)*(n-2),1);        % Se inicializa rhs con ceros.
phi_approx = zeros(n,m);           % Se inicializa phi_approx

% Ensamblamos la Matriz A
dB = diag(4*ones(1,(n-2)));        % Hacemos una matriz diagonal.
dBp1 = diag(1*ones(1,(n-2)-1),1);  % Creamos la matriz diagonal superior.
dBm1 = diag(1*ones(1,(n-2)-1),-1); % Creamos la matriz diagonal inferior.
B = (dB - dBp1 - dBm1);            % Juntamos las matrices para obtener B.
I = -eye(n-2);                     % Creamos la matriz identidad.
temp = 1;                          % Creamos una variable temporal.
for i=1:(n-2):(m-2)*(n-2)          % Se recorre de 1 a (m-2)*(m-2) dando salida
    A(i:temp*(n-2),i:temp*(n-2)) = B;% Se asignan los bloques con B.
    if temp*(n-2) < (m-2)*(n-2)    % Se buscan los espacios para los bloques
        A(temp*(n-2)+1:temp*(n-2)+(n-2),i:temp*(n-2)) = I;% Se agrega el bloque de I a la izquierda
        A(i:temp*(n-2),temp*(n-2)+1:temp*(n-2)+(n-2)) = I;% Se agrega el bloque de I a la derecha.
    end
    temp = temp + 1;               % Se hace un incremento en temp.
end
%disp(A)                                 % Imprimimos A en pantalla.
%spy(A)
%disp(x)
%disp(x(1,2))
%disp(x(2,1))

% Agregamos condiciones de frontera
for i =2:n-1
    temp = (i-1);                  % Se cambia el valor de temp.
    rhs(temp) = rhs(temp) + phi(x(i,1),y(i,1));% Se agrega la condición inicial al
    temp = (i-1) + (n-2)*(m-3);% Se cambia el valor de temp.
    %disp(temp)
    rhs(temp) = rhs(temp) + phi(x(i,m),y(i,m));% Se agrega la condición inicial al
end
for j = 2:m-1
    temp = (j-2)*(n-2)+1;    % Se cambia el valor de temp.
    rhs(temp) = rhs(temp) + phi(x(1,j),y(1,j));% Se agrega la condicion inicial al
    temp = (j-1)*(n-2);% Se cambia el valor de temp.
    %disp(temp)
    rhs(temp) = rhs(temp) + phi(x(n,j),y(n,j));% Se agrega la condicion incial al 
end
for i=2:m-1
    for j = 2:n-1
        temp = i+j-3+(n-3)*(i-2);% Se cambia el valor de temp.
        %disp(temp)
        rhs(temp) = rhs(temp) - h^2*f(x(i,j),y(i,j));% Se agrega la condicion inicial al
    end
end

%%%%%%%%%% Características de la matriz %%%%%%%%%%%
%C = cond(A);
%fprintf('El número de condición es: ');
%disp(C)

%addpath('C:\Users\Dell\Documents\Curso_Computo_Cientifico\Poisson\examples\caracteristicas_Matrices');
%isdom = IsDiagDom(A);
%if isdom == 0
%    disp ('Matrix A is not diagonally-dominant');
%elseif isdom == 1
%        disp ('Matrix A is diagonally-dominant'); 
%end

%%%%% Resolvemos el sistema lineal %%%%%%%%%%%%%%
%tStart = cputime;
%u = A\rhs;                        % Utiliza mldivide para resolver el sistema
%tiempo = cputime - tStart;

tStart = cputime;                 %% Inicia método LU %%%%%%%
[L,U,P] = lu(A);                  % Factorizamos la matriz en una triangular inferior y una superior
v = L\(P*rhs);                    % Resolvemos la matriz inferior
u = U\v;                          % Resolvemos la matriz superior a partir del resultado anterior
tiempo = cputime - tStart;

%tStart = cputime;                 %% Inicia método Cholesky %%%%%%%
%R = chol(A);                      % Factorizamos la matriz en una triangular inferior
%u = R\(R'\rhs);                   % Resolvemos 
%tiempo = cputime - tStart;
%disp(u)

%addpath('C:\Users\Dell\Documents\Curso_Computo_Cientifico\Poisson\examples\metodos_iterativos');
%tStart = cputime;
%x0 = zeros((m-2)*(n-2),1);        % Vector inicial.
%u = GaussSeidel(A,rhs,x0,(m-2)*(n-2),0.000001,500);  %% Método implementado por mí
%u = Gauss_Siedel(A,rhs,x0,0.000001);   %% Implementación descargada de internet(https://www.mathworks.com/matlabcentral/fileexchange/73488-gauss-seidel-iterative-method)
%tiempo = cputime - tStart;

%tStart = cputime;
%x0 = zeros((m-2)*(n-2),1);        % Vector inicial.
%u = Jacobi(A,rhs,x0,(m-2)*(n-2),0.000001,500);
%tiempo = cputime - tStart;

%disp(u) 
%% Guardamos la solución
utemp = reshape(u,n-2,m-2);        % Cambiamos de vector a matriz
for i = 1:n
    phi_approx(i,1) = phi(x(i,1),y(i,1));% Se agrega la condicion de frontera
    phi_approx(i,m) = phi(x(i,m),y(i,m));% Se agrega la condicion de frontera
end
%disp(phi_approx(:,1))
%disp(phi_approx(:,m))
for j = 2:m-1
    phi_approx(1,j) = phi(x(1,j),y(1,j));% Se agrega la condicion de frontera
    %disp(phi_approx(1,j))
    %disp(j)
    phi_approx(n,j) = phi(x(n,j),y(n,j));% Se agrega la condicion de frontera
end
%disp(phi_approx(1,:))
%disp(phi_approx(n,:))
phi_approx(2:(n-1),2:(m-1)) = utemp;       % Se completa la matriz de phi_approx.

%disp(phi_approx)
%disp(x)
%disp(y)
% Calculamos la solucion exacta
phi_exacta = phi(x,y);
%disp(phi_exacta)

%% Graficamos la solución
% scrsz = get(groot,'ScreenSize'); % Se obtienen los limites de la
% figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]);% Se crea la figura del
% tamaño de
%figure                             
%subplot(1,2,1)                     % Se divide la grafica en 2.
%surf(x,y,phi_approx);              % Se grafica la solucion aproximada.
%title('Aproximacion')              % Se agrega titulo a la grafica
%subplot(1,2,2)                     % Se usa la otra parte de la grafica.
%surf(x,y,phi_exacta);              % Se grafica la solucion exacta.
%title('Solucion exacta')           % Se agrega titulo a la grafica
end