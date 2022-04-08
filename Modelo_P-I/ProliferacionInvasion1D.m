function [c_approx,x,cont11] = ProliferacionInvasion1D(n,t,g,v)
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
%   n     Entero       Numero de nodos en la malla en el espacio en x.
%   t     Entero       Numero de pasos en el tiempo
%   v     Real         Coeficiente de difusion
%   f     funcion      funcion para condicion inicial y de frontera
%
% Argumentos de salida
%   c_approx     matriz       Matriz con la solucion aproximada
%   x            matriz       Vector de coordenadas x de la malla
%   t            real         Tiempo total (dias, minutos, segundos)
% 
% Ejemplo de uso:
%       [c_approx,x,tiempoEjec,iter] = ProliferacionInvasion2D(30,5,@g,0.0013);
%
%
%% Se inicializan las variables
close all                            % Cerramos todas las ventanas de figuras.
x = linspace(0,n,n);                 % Se crea la malla espacial en x.
dt = 1;                            % Tamaño de paso de integracion en el tiempo.
dx = x(2) - x(1);                % Se calcula dx.
c_approx = zeros(n);                % Se inicializa c_approx con ceros. Aquí se guarda las aproximaciones.
%% P A R A    A R M A R    M A T R I C E S %%%%%%%%
A = zeros((n-2),(n-2));% Se inicializa A con ceros.
rhs = zeros((n-2),1);        % Se inicializa rhs con ceros.
%stb = abs(a*dt/dx + b*dt/dy);        % Calculamos c para ver la estabilidad.
%if stb > 1/2
%    fprintf('El codigo puede no funcionar c = %f\n',c);
%end
%% P A R A   A P L I C A R   M É T O D O   I T E R A T I V O  %%%%%
%err = 1;                              % Se incializa err con 1.
%tol = sqrt(eps);                      % Se impone una tolerancia para err.
%tol = 0.000001;

%% Agregamos condiciones inicial para el tiempo T=1
for i = 1:n
        c_approx(i) = 4000;
end

%% Agregamos la condicion de frontera.  (Para el tiempo t=1)
c_approx(1) = 0;
c_approx(n) = 0;


%% Llenamos la matriz A para llevar nuestra discretización a la forma Ax=b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   E S T O      P A R A       R E S O L V E R      C O N     %%%%%%%
%%%      M É T O D O S      D I R E C T O S                     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensamblamos la Matriz A
ax = v/(dx)^2;
by = v/(dy)^2;
dB = diag(-(2*ax+2*by)*ones(1,(n-2)));        % Hacemos una matriz diagonal.
dBp1 = diag(by*ones(1,(n-2)-1),1);  % Creamos la matriz diagonal superior.
dBm1 = diag(by*ones(1,(n-2)-1),-1); % Creamos la matriz diagonal inferior.
B = (dB + dBp1 + dBm1);            % Juntamos las matrices para obtener B.
Ix = ax*eye(n-2);                     % Creamos la matriz identidad.
temp = 1;                          % Creamos una variable temporal.
for i=1:(n-2):(m-2)*(n-2)          % Se recorre de 1 a (m-2)*(m-2) dando salida
    A(i:temp*(n-2),i:temp*(n-2)) = B;% Se asignan los bloques con B.
    if temp*(n-2) < (m-2)*(n-2)    % Se buscan los espacios para los bloques
        A(temp*(n-2)+1:temp*(n-2)+(n-2),i:temp*(n-2)) = Ix;% Se agrega el bloque de I a la izquierda
        A(i:temp*(n-2),temp*(n-2)+1:temp*(n-2)+(n-2)) = Ix;% Se agrega el bloque de I a la derecha.
    end
    temp = temp + 1;               % Se hace un incremento en temp.
end


%% B O R R A R    S I     E S      N E C E S A R I O  %%%%%%%%%%%%%% 
%temp = 1;


%for i = 1:(n-2)*(m-2):(m-2)*(n-2)*(t-2)
%    A(i:temp*(m-2)*(n-2),i:temp*(n-2)*(m-2)) = Ai;
%    if temp*(n-2)*(m-2) < (m-2)*(n-2)*(t-2)    % Se buscan los espacios para los bloques
%        A(temp*(n-2)*(m-2)+1:temp*(n-2)*(m-2)+(n-2)*(m-2),i:temp*(n-2)*(m-2)) = It;% Se agrega el bloque de I a la izquierda
%        A(i:temp*(n-2)*(m-2),temp*(n-2)*(m-2)+1:temp*(n-2)*(m-2)+(n-2)*(m-2)) = It;% Se agrega el bloque de I a la derecha.
%    end
%    temp = temp + 1;
%end
%disp(A)                                 % Imprimimos A en pantalla.
%spy(A)


%% Hacemos el vector b del sistema de ecuaciones algebraicas

for k = 1:dt:t                         %ciclo para integración de tiempo
    
    % Armar rhs
    
    for i =2:n-1
        temp = (i-1);                              % Se cambia el valor de temp.
        rhs(temp) = rhs(temp) + by*c_approx(1,i);  % Se agrega la condición inicial al vertor b
        temp = (i-1) + (n-2)*(m-3);                % Se cambia el valor de temp.
        %disp(temp)
        rhs(temp) = rhs(temp) + by*c_approx(m,i);  % Se agrega la condición inicial al vector b
    end
    for j = 2:m-1
        temp = (j-2)*(n-2)+1;                      % Se cambia el valor de temp.
        rhs(temp) = rhs(temp) + ax*c_approx(j,1);  % Se agrega la condicion inicial al vector b
        temp = (j-1)*(n-2);                        % Se cambia el valor de temp.
       %disp(temp)
        rhs(temp) = rhs(temp) + ax*c_approx(j,n);  % Se agrega la condicion incial al vector b
    end
    for i=2:m-1
        for j = 2:n-1
            temp = i+j-3+(n-3)*(i-2);               % Se cambia el valor de temp.
           %disp(temp)
            rhs(temp) = rhs(temp) - g(x(i,j),y(i,j),c_approx(i,j));% Se agrega la condicion inicial al vector b
        end
    end
    
%% Resolvemos el sistema algebraico %%
    %tStart = cputime;                 %% Inicia método LU %%%%%%%
    [L,U,P] = lu(A);                  % Factorizamos la matriz en una triangular inferior y una superior
    v = L\(P*rhs);                    % Resolvemos la matriz inferior
    u = U\v;                          % Resolvemos la matriz superior a partir del resultado anterior
    %tiempo = cputime - tStart;

%% Guardamos la solución        
    utemp = reshape(u,n-2,m-2);        % Cambiamos de vector a matriz
    c_approx(2:(n-1),2:(m-1)) = utemp;       % Se completa la matriz de phi_approx.


%% Imprimir solucion
    plot3(x,y,c_approx);

end
xlabel('x')
ylabel('y')
zlabel('c_approx')
title('Aproximacion')

%disp(c_approx)
%% Hacemos el metodo para Diferencias de segundo orden   (Esta solución es un tipo Gauss-Seidel sin armar una matriz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% E S T O   P A R A    R E S O L V E R L O    CON   %%%%%
%%%% U N    M É T O D O    I T E R A T I V O           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont11 = 1;
%tStart = cputime;
%while (err >= tol) && (cont11 < 1500)
%while (err >= tol) 
%    err=0;
%    for k = 2:t
%        for i = 2:m-1
%            for j = 2:n-1
               %act = dt*(g(x(i,j),y(i,j),c_approx(i,j,k-1),k-1) + (v/dx^2)*(c_approx(i+1,j,k-1) - 2*c_approx(i,j,k-1) + c_approx(i-1,j,k-1)) + (v/dy^2)*(c_approx(i,j+1,k-1) - 2*c_approx(i,j,k-1) + c_approx(i,j-1,k-1))) - c_approx(i,j,k-1); % Se hace el llenado de u para todos los puntos.
%                act = (1/((2*v/dx^2)+(2*v/dy^2)-(1/dt)))*((v/dx^2)*(c_approx(i+1,j,k-1)+c_approx(i-1,j,k-1))+ (v/dy^2)*(c_approx(i,j+1,k-1)+c_approx(i,j-1,k-1)) + g(x(i,j),y(i,j),c_approx(i,j,k-1),k-1));
%                err = max(err,abs(act-c_approx(i,j,k-1)));% Se calcula un nuevo err.
%                c_approx(i,j,k) = act;
%            end
%        end
%    end
%    cont11=cont11+1;
%    disp(cont11)
%    disp(err)
%end
%tiempo = cputime - tStart;

%% Graficamos la solucion aproximada   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% E S T A    P A R T E    N O    H A     S I D O    %%%%
%%% R E V I S A D A    P O R    E L     A S E S O R   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
