function [phi_approx, phi_exacta,x,y,z,tiempo] = Poisson3D(m,n,p,phi,f)
% Función que calcula una aproximación a la solución de la ecuación de
% Possion en 2D.
%
%phi_xx + phi_yy + phi_zz = -fhjvg
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
% [phi_aprox,phi_exacta,x,y,z,tiempo] = Poisson3D(11,11,11,@phi,@f);
%
% Inicializamos las variables
close all                          % Cierra ventanas de figuras abiertas
x = linspace(0,1,m);              % Se crea la discretización en x
y = linspace(0,1,n);              % Se crea la discretización en y
z = linspace(0,1,p);              % Se crea la discretización en z
h = x(2) - x(1);                   % Se calcula h
[x,y,z] = meshgrid(x,y,z);             % Creamos la malla completa.
Ai = zeros((m-2)*(n-2),(m-2)*(n-2));% Se inicializa A con ceros.
A = zeros((m-2)*(n-2)*(p-2),(m-2)*(n-2)*(p-2));% Se inicializa A con ceros.
rhs = zeros((m-2)*(n-2)*(p-2),1);        % Se inicializa rhs con ceros.
phi_approx = zeros(n,m,p);           % Se inicializa phi_approx

% Ensamblamos la Matriz A
dB = diag(4*ones(1,(n-2)));        % Hacemos una matriz diagonal.
dBp1 = diag(1*ones(1,(n-2)-1),1);  % Creamos la matriz diagonal superior.
dBm1 = diag(1*ones(1,(n-2)-1),-1); % Creamos la matriz diagonal inferior.
B = (dB - dBp1 - dBm1);            % Juntamos las matrices para obtener B.
I = -eye(n-2);                     % Creamos la matriz identidad.
temp = 1;                          % Creamos una variable temporal.
for i=1:(n-2):(m-2)*(n-2)          % Se recorre de 1 a (m-2)*(m-2) dando salida
    Ai(i:temp*(n-2),i:temp*(n-2)) = B;% Se asignan los bloques con B.
    if temp*(n-2) < (m-2)*(n-2)    % Se buscan los espacios para los bloques
        Ai(temp*(n-2)+1:temp*(n-2)+(n-2),i:temp*(n-2)) = I;% Se agrega el bloque de I a la izquierda
        Ai(i:temp*(n-2),temp*(n-2)+1:temp*(n-2)+(n-2)) = I;% Se agrega el bloque de I a la derecha.
    end
    temp = temp + 1;               % Se hace un incremento en temp.
end
Ip = -eye((n-2)*(m-2));                     % Creamos la matriz identidad.
temp = 1;
for i = 1:(n-2)*(m-2):(m-2)*(n-2)*(p-2)
    A(i:temp*(m-2)*(n-2),i:temp*(n-2)*(m-2)) = Ai;
    if temp*(n-2)*(m-2) < (m-2)*(n-2)*(p-2)    % Se buscan los espacios para los bloques
        A(temp*(n-2)*(m-2)+1:temp*(n-2)*(m-2)+(n-2)*(m-2),i:temp*(n-2)*(m-2)) = Ip;% Se agrega el bloque de I a la izquierda
        A(i:temp*(n-2)*(m-2),temp*(n-2)*(m-2)+1:temp*(n-2)*(m-2)+(n-2)*(m-2)) = Ip;% Se agrega el bloque de I a la derecha.
    end
    temp = temp + 1;
end
%disp(A)                                 % Imprimimos A en pantalla.
%d = eig(A);
%disp(d)
%disp(x(1,2))
%disp(x(2,1))

% Agregamos condiciones de frontera

for j =2:n-1
    for k = 2:p-1
        temp1 = (k-1);                  % Se cambia el valor de temp.
        rhs(temp1) = rhs(temp1) + phi(x(1,j,k),y(1,j,k),z(1,j,k));% Se agrega la condición inicial al
        temp2 = (k-1) + (n-2)*(m-2)*(p-2) - (n-2)*(p-2);% Se cambia el valor de temp.
        rhs(temp2) = rhs(temp2) + phi(x(m,j,k),y(m,j,k),z(m,j,k));% Se agrega la condición inicial al
    end
end
for i = 2:m-1
    for k=2:p-1
        temp1 = (i-2)*(n-2)*(p-2) + (k-1);    % Se cambia el valor de temp.
        rhs(temp1) = rhs(temp1) + phi(x(i,1,k),y(i,1,k),z(i,1,k));% Se agrega la condicion inicial al
        temp2 = (i-1)*(n-2)*(p-2) - (p-2) + (k-1);% Se cambia el valor de temp.
        rhs(temp2) = rhs(temp2) + phi(x(i,n,k),y(i,n,k),z(i,n,k));% Se agrega la condicion incial al 
    end
end
for i=2:m-1
    for j=2:n-1
        temp1 = (i-2)*(m-2)*(n-2) + (n-2)*(j-2) + 1;    % Se cambia el valor de temp.
        rhs(temp1) = rhs(temp1) + phi(x(i,j,1),y(i,j,1),z(i,j,1));% Se agrega la condicion inicial al
        temp2 = (i-2)*(m-2)*(n-2) + (n-2)*(j-1);% Se cambia el valor de temp.
        rhs(temp2) = rhs(temp2) + phi(x(i,j,p),y(i,j,p),z(i,j,p));% Se agrega la condicion incial al 
    end
end
for i=2:m-1
    for j = 2:n-1
        for k = 2:p-1
            temp = i+j+k-5+(n-3)*(j-2)+((n-2)*(m-2)-1)*(i-2);% Se cambia el valor de temp.
            %disp(temp)
            rhs(temp) = rhs(temp) - h^2*f(x(i,j,k),y(i,j,k),z(i,j,k));% Se agrega la condicion inicial al
        end
    end
end

%%%%% Resolvemos el sistema lineal %%%%%%%%%%%%%%
tStart = cputime;
u = A\rhs;                        % Utiliza mldivide para resolver el sistema
tiempo = cputime - tStart;

%tStart = cputime;                 %% Inicia método LU %%%%%%%
%[L,U,P] = lu(A);                  % Factorizamos la matriz en una triangular inferior y una superior
%v = L\(P*rhs);                    % Resolvemos la matriz inferior
%u = U\v;                          % Resolvemos la matriz superior a partir del resultado anterior
%tiempo = cputime - tStart;

%tStart = cputime;                 %% Inicia método Cholesky %%%%%%%
%R = chol(A);                      % Factorizamos la matriz en una triangular inferior
%u = R\(R'\rhs);                   % Resolvemos 
%tiempo = cputime - tStart;
%disp(u)
% Guardamos la solución
utemp = reshape(u,m-2,n-2,p-2);        % Cambiamos de vector a matriz
for i = 1:m
    for j = 1:n
        phi_approx(i,j,1) = phi(x(i,j,1),y(i,j,1),z(i,j,1));% Se agrega la condicion de frontera
        phi_approx(i,j,p) = phi(x(i,j,p),y(i,j,p),z(i,j,p));% Se agrega la condicion de frontera
    end
end
%disp(phi_approx(:,1))
%disp(phi_approx(:,m))
for k = 2:p-1
    phi_approx(1,1,k) = phi(x(1,1,k),y(1,1,k),z(1,1,k));% Se agrega la condicion de frontera
    phi_approx(1,n,k) = phi(x(1,n,k),y(1,n,k),z(1,n,k));
    phi_approx(m,1,k) = phi(x(m,1,k),y(m,1,k),z(m,1,k));
    phi_approx(m,n,k) = phi(x(m,n,k),y(m,n,k),z(m,n,k));
end
%disp(phi_approx(1,:))
%disp(phi_approx(n,:))
phi_approx(2:(m-1),2:(n-1),2:(p-1)) = utemp;       % Se completa la matriz de phi_approx.

%disp(phi_approx)
%disp(x)
%disp(y)
% Calculamos la solucion exacta
phi_exacta = phi(x,y,z);
%disp(phi_exacta)

% Graficamos la solución
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