function [c_approx,x,y,tiempo] = ProliferacionInvasion2D(m,n,t,f,g,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretizacion mediante el método de diferencias finitas %%
% Para el modelo de proliferacion invasion.                %%
% Sandra I. García Mendoza                                 %%
% Noviembre 2021                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  c_t = v*(c_xx + c_yy) + g   % f = c(X,t): representa la densidad del tumor
%                              % g: representa la función que modela la
%                              invasión y el Tratamiento en caso de tener alguno.
%
% Se hace un metodo iterativo para la aproximacion.
%
% Argumentos de entrada
%   m     Entero       Numero de nodos en la malla en el espacio en x.
%   n     Entero       Numero de nodos en la malla en el espacio en y.
%   t     Entero       Numero de pasos en el tiempo
%   v     Real         Coeficiente de difusion
%   f     funcion      funcion para condicion inicial y de frontera
%
% Argumentos de salida
%   c_approx     matriz       Matriz con la solucion aproximada
%   x            matriz       Vector de coordenadas x de la malla
%   y            matriz       Vector de coordenadas y de la malla
% 
% Ejemplo de uso:
%       [c,x,y] = ProliferacionInvasion2D(21,21,10,@f,@g,0.2);
%
%
% Se inicializan las variables
close all                            % Cerramos todas las ventanas de figuras.
T = linspace(0,90,t);                 % Se crea la malla temporal.
x = linspace(0,200,m);                 % Se crea la malla espacial en x.
y = linspace(0,200,n);                 % Se crea la malla espacial en y.
[x,y] = meshgrid(x,y);               % Se crea la malla espacial.
dt = T(2) - T(1);                    % Se calcula dt.
dx = x(1,2) - x(1,1);                % Se calcula dx.
dy = y(2,1) - y(1,1);                % Se calcula dx.
%c_exacta = zeros(m,n,t);             % Se inicializa u_exacta con ceros.
c_approx = zeros(m,n,t);             % Se inicializa u_approx con ceros.
A = zeros((m-2)*(n-2)*(t-2),(m-2)*(n-2)*(t-2));% Se inicializa A con ceros.
Ai = zeros((m-2)*(n-2),(m-2)*(n-2));% Se inicializa A con ceros.
rhs = zeros((m-2)*(n-2)*(t-2),1);        % Se inicializa rhs con ceros.
%stb = abs(a*dt/dx + b*dt/dy);        % Calculamos c para ver la estabilidad.
%if stb > 1/2
%    fprintf('El codigo puede no funcionar c = %f\n',c);
%end

%% Agregamos condiciones inicial   para el tiempo T=1
for i = 1:m
    for j = 1:n
        c_approx(i,j,1) = f(x(i,j),y(i,j),T(1));  % Se agrega la condicion inicial.
    end
end

%% Agregamos la condicion de frontera.  (Para los tiempos t>=2)
for k = 2:t
    for i = 1:n
        c_approx(i,1,k) = f(x(i,1),y(i,1),T(k));  % Se agrega la condicion de frontera en x=0.
        c_approx(i,n,k) = f(x(i,m),y(i,m),T(k));  % Se agrega la condicion de frontera en x=m.
    end
    for j = 2:m
        c_approx(1,j,k) = f(x(1,j),y(1,j),T(k));  % Se agrega la condicion de frontera en y=0.
        c_approx(m,j,k) = f(x(n,j),y(n,j),T(k));  % Se agrega la condicion de frontera en y=n.
    end
end

%% Llenamos la matriz A para llevar nuestra discretización a la forma Ax=b

% Ensamblamos la Matriz A
ax = v/(dx)^2;
by = v/(dy)^2;
ct = 1/dt;
dB = diag((2*ax+2*by-ct)*ones(1,(n-2)));        % Hacemos una matriz diagonal.
dBp1 = diag(-by*ones(1,(n-2)-1),1);  % Creamos la matriz diagonal superior.
dBm1 = diag(-by*ones(1,(n-2)-1),-1); % Creamos la matriz diagonal inferior.
B = (dB - dBp1 - dBm1);            % Juntamos las matrices para obtener B.
Ix = -ax*eye(n-2);                     % Creamos la matriz identidad.
It = ct*eye((n-2)*(m-2));              % Creamos la matriz identidad.
temp = 1;                          % Creamos una variable temporal.
for i=1:(n-2):(m-2)*(n-2)          % Se recorre de 1 a (m-2)*(m-2) dando salida
    Ai(i:temp*(n-2),i:temp*(n-2)) = B;% Se asignan los bloques con B.
    if temp*(n-2) < (m-2)*(n-2)    % Se buscan los espacios para los bloques
        Ai(temp*(n-2)+1:temp*(n-2)+(n-2),i:temp*(n-2)) = Ix;% Se agrega el bloque de I a la izquierda
        Ai(i:temp*(n-2),temp*(n-2)+1:temp*(n-2)+(n-2)) = Ix;% Se agrega el bloque de I a la derecha.
    end
    temp = temp + 1;               % Se hace un incremento en temp.
end
temp = 1;
for i = 1:(n-2)*(m-2):(m-2)*(n-2)*(t-2)
    A(i:temp*(m-2)*(n-2),i:temp*(n-2)*(m-2)) = Ai;
    if temp*(n-2)*(m-2) < (m-2)*(n-2)*(t-2)    % Se buscan los espacios para los bloques
        A(temp*(n-2)*(m-2)+1:temp*(n-2)*(m-2)+(n-2)*(m-2),i:temp*(n-2)*(m-2)) = It;% Se agrega el bloque de I a la izquierda
        A(i:temp*(n-2)*(m-2),temp*(n-2)*(m-2)+1:temp*(n-2)*(m-2)+(n-2)*(m-2)) = It;% Se agrega el bloque de I a la derecha.
    end
    temp = temp + 1;
end
%disp(A)                                 % Imprimimos A en pantalla.
%spy(A)


% Hacemos el vector b del sistema de ecuaciones algebraicas

for k = 2:t-2                         %c_approx(m,n,t)ijk
    for i =2:n-1
        temp = (i-1);                              % Se cambia el valor de temp.
        rhs(temp) = rhs(temp) + by*c_approx(1,i,t);%phi(x(i,1),y(i,1));% Se agrega la condición inicial al vertor b
        temp = (i-1) + (n-2)*(m-3);                % Se cambia el valor de temp.
        %disp(temp)
        rhs(temp) = rhs(temp) + by*c_approx(m,i,t);%phi(x(i,m),y(i,m));% Se agrega la condición inicial al vector b
    end
    for j = 2:m-1
        temp = (j-2)*(n-2)+1;                      % Se cambia el valor de temp.
        rhs(temp) = rhs(temp) + ax*c_approx(j,1,t);%phi(x(1,j),y(1,j));% Se agrega la condicion inicial al vector b
        temp = (j-1)*(n-2);                        % Se cambia el valor de temp.
        %disp(temp)
        rhs(temp) = rhs(temp) + ax*c_approx(j,n,t);%phi(x(n,j),y(n,j));% Se agrega la condicion incial al vector b
    end
    for i=2:m-1
        for j = 2:n-1
            temp = i+j-3+(n-3)*(i-2);               % Se cambia el valor de temp.
            %disp(temp)
            rhs(temp) = rhs(temp) - g(x(i,j),y(i,j),@f,t);% Se agrega la condicion inicial al vector b
        end
    end
end

%% Resolvemos el sistema algebraico %%
tStart = cputime;                 %% Inicia método LU %%%%%%%
[L,U,P] = lu(A);                  % Factorizamos la matriz en una triangular inferior y una superior
v = L\(P*rhs);                    % Resolvemos la matriz inferior
u = U\v;                          % Resolvemos la matriz superior a partir del resultado anterior
tiempo = cputime - tStart;

%% Guardamos la solución        %%% R E V I S A R   C O M O    Q U E D A   P A R A   E S T E   C A S O %%%%%
utemp = reshape(u,n-2,m-2,t-2);        % Cambiamos de vector a matriz
%for i = 1:n
%    phi_approx(i,1) = phi(x(i,1),y(i,1));% Se agrega la condicion de frontera
%    phi_approx(i,m) = phi(x(i,m),y(i,m));% Se agrega la condicion de frontera
%end
%disp(phi_approx(:,1))
%disp(phi_approx(:,m))
%for j = 2:m-1
%    phi_approx(1,j) = phi(x(1,j),y(1,j));% Se agrega la condicion de frontera
    %disp(phi_approx(1,j))
    %disp(j)
%    phi_approx(n,j) = phi(x(n,j),y(n,j));% Se agrega la condicion de frontera
%end
%disp(phi_approx(1,:))
%disp(phi_approx(n,:))


c_approx(2:(n-1),2:(m-1),2:(t-1)) = utemp;       % Se completa la matriz de phi_approx.


%disp(c_approx)
%% Hacemos el metodo para Diferencias de segundo orden   (Esta solución es un tipo Gauss-Seidel)
%tStart = cputime;
%for k = 2:t
%    for i = 2:m-1
%        for j = 2:n-1
%            c_approx(i,j,k) = c_approx(i,j,k-1) + dt*((v/dx^2)*(c_approx(i+1,j,k-1) - 2*c_approx(i,j,k-1) + c_approx(i-1,j,k-1)) + (v/dy^2)*(c_approx(i,j+1,k-1) - 2*c_approx(i,j,k-1) + c_approx(i,j-1,k-1))); % Se hace el llenado de u para todos los puntos.
%        end
%    end
%end
%tiempo = cputime - tStart;
%% Calculamos la solucion exacta  %%%% Escribir la ecuacion completa en la forma inicial
%for k = 1:t
%    for i = 1:m
%        for j = 1:n
%            c_exacta(i,j,k) = f(x(i,j),y(i,j),T(k));
%        end
%    end
%end

%% Graficamos la solucion 
%scrsz = get(groot, 'ScreenSize');  % Se obtienen los limites de la
%figure('OuterPosition',[1 1 scrsz(3) scrsz(4)]); % Se crea una figura del tamaño
%figure
%p = ceil(t/100);                     % Se calcula p para no graficar todos los 
%for i = 1:p:t-p
%    subplot(2,5,1)
%    surf(x,y,c_approx(:,:,1))                   % Se grafica la aproximación de c.
%    title('Aproximacion')
    %axis([0 1 0 1 -1 1])                % Se ajustan los limites de la grafica.
    %pause(0.8)
%    subplot(2,5,2)
%    surf(x,y,c_approx(:,:,2))                   % Se grafica la solución exacta c.
%    subplot(2,5,3)
%    surf(x,y,c_approx(:,:,3))  
%    subplot(2,5,4)
%    surf(x,y,c_approx(:,:,4))  
%    subplot(2,5,5)
%    surf(x,y,c_approx(:,:,5))  
%    subplot(2,5,6)
%    surf(x,y,c_approx(:,:,6))  
%    subplot(2,5,7)
%    surf(x,y,c_approx(:,:,7))  
%    subplot(2,5,8)
%    surf(x,y,c_approx(:,:,8))  
%    subplot(2,5,9)
%    surf(x,y,c_approx(:,:,9))  
%    subplot(2,5,10)
%    surf(x,y,c_approx(:,:,10))  
    %title('Solucion exacta')
    %axis([0 1 0 1 -1 1])                % Se ajustan los limites de la grafica.
    %pause(0.01)
%end
end
