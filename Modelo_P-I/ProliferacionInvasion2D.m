function [c_approx,x,y,cont11] = ProliferacionInvasion2D(m,n,t,g,v)
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
%   t            real         Tiempo total (dias, minutos, segundos)
% 
% Ejemplo de uso:
%       [c_approx,x,y,tiempoEjec,iter] = ProliferacionInvasion2D(10,10,5,@g,0.0.0013);
%
%
%% Se inicializan las variables
close all                            % Cerramos todas las ventanas de figuras.
x = linspace(0,100,m);                 % Se crea la malla espacial en x.
y = linspace(0,100,n);                 % Se crea la malla espacial en y.
[x,y] = meshgrid(x,y);               % Se crea la malla espacial.
dt = 0.001;                            % Tamaño de paso de integracion en el tiempo.
dx = x(1,2) - x(1,1);                % Se calcula dx.
dy = y(2,1) - y(1,1);                % Se calcula dx.
c_approx = zeros(m*n,1);                % Se inicializa c_approx con ceros. Aquí se guarda las aproximaciones.
%% P A R A    A R M A R    M A T R I C E S %%%%%%%%
A = zeros(m*n,m*n);% Se inicializa A con ceros.
rhs = zeros(m*n,1);        % Se inicializa rhs con ceros.
stb = abs(v*dt/dx^2);        % Calculamos c para ver la estabilidad.
if stb > 1/2
    fprintf('El codigo puede no funcionar c = %f\n',c);
end
%% P A R A   A P L I C A R   M É T O D O   I T E R A T I V O  %%%%%
%err = 1;                              % Se incializa err con 1.
%tol = sqrt(eps);                      % Se impone una tolerancia para err.
%tol = 0.000001;

%%  % Calculamos c para ver la estabilidad.
cest = abs(v*dt/dx^2);                   
if cest > 1/2
    fprintf('El codigo puede no funcionar c = %f\n',cest);
end

%% Agregamos condiciones inicial para el tiempo T=1
for i=1:m-1                        % condiciones normales i=1:m-2  %condiciones auxiliares i = m/2:(m/2)+2 
    for j=2:n-1                     % j=2:n-1                       %j = n/2:(n/2)+2
        c_approx((m*i)+j) = 4000;
    end
end
%c_approx = c_approx + 4000;


%% Llenamos la matriz A para llevar nuestra discretización a la forma Ax=b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   E S T O      P A R A       R E S O L V E R      C O N     %%%%%%%
%%%      M É T O D O S      D I R E C T O S                     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensamblamos la Matriz A
ax = v/(dx)^2;
by = v/(dy)^2;
dB = diag((2*ax+2*by)*ones(1,(n)));        % Hacemos una matriz diagonal.
dBp1 = diag(-by*ones(1,(n)-1),1);  % Creamos la matriz diagonal superior.
dBm1 = diag(-by*ones(1,(n)-1),-1); % Creamos la matriz diagonal inferior.
B = (dB - dBp1 - dBm1);            % Juntamos las matrices para obtener B.
Ix = ax*eye(n);                     % Creamos la matriz identidad.
temp = 1;                          % Creamos una variable temporal.
for i=1:(n):(m)*(n)          % Se recorre de 1 a (m-2)*(m-2) dando salida
    A(i:temp*(n),i:temp*(n)) = B;% Se asignan los bloques con B.
    if temp*(n) < (m)*(n)    % Se buscan los espacios para los bloques
        A(temp*(n)+1:temp*(n)+(n),i:temp*(n)) = Ix;% Se agrega el bloque de I a la izquierda
        A(i:temp*n,temp*n+1:temp*n+n) = Ix;% Se agrega el bloque de I a la derecha.
    end
    temp = temp + 1;               % Se hace un incremento en temp.
end

figure
u_approx = reshape(c_approx,n,m);
disp(u_approx)
disp(c_approx)
subplot(1,3,1)
surf(x,y,u_approx);
xlabel('x')
ylabel('y')
zlabel('c_approx')
%hold on 

%disp(A)                                 % Imprimimos A en pantalla.
%spy(A)



%% Actualizacion de Condiciones de frontera
%for i = 1:m
%        if (i ~= 1) || (i ~= m)
%            c_approx(m*(i-1)+1) = 0;
%            c_approx(m*(i-1)+n) = 0;
%        elseif i == 1
%            for j = 1:n
%                c_approx(j) = 0;
%            end
%        elseif i == m
%            for j = 1:n
%                c_approx((m*n - n)+j) = 0;
%            end
%        end
%end
%% Paso de integración de tiempo

for k = 0:dt:t                       

 % Ensamblado del lado derecho ( vector rhs )
    
    %for i=1:n*m
    %        rhs(i) = v*c_approx(i) - g(x,y,c_approx(i)) ;% Se agrega la condicion inicial al vector b
    %end   
    rhs = v*c_approx;
   
    %disp(rhs)
    
%% Resolvemos el sistema algebraico %%
    %tStart = cputime;                 %% Inicia método LU %%%%%%%
    [L,U,P] = lu(A);                  % Factorizamos la matriz en una triangular inferior y una superior
    vM = L\(P*rhs);                    % Resolvemos la matriz inferior
    du = U\vM;                          % Resolvemos la matriz superior a partir del resultado anterior
    %tiempo = cputime - tStart;
    %disp(du)


%% Integración de tiempo Método de Euler
    c_approx = c_approx + dt*du;

    
%% Integración de tiempo Runge Kutta 2
    % k1
%    wi_approx = c_approx + (dt/2)*du;
    
    % Resolver nuevamente con la evaluación en t+h/2,wi_approx
    % Actualizamos rhs (vector b)
%    for i=1:n*m
%            rhs(i) = - g(x,y,wi_approx(i)) ;% Se agrega la condicion inicial al vector b
%    end
    % Resolvemos 
%    [L,U,P] = lu(A);                  % Factorizamos la matriz en una triangular inferior y una superior
%    vi = L\(P*rhs);                    % Resolvemos la matriz inferior
%    duh2 = U\vi;   
    
    % k2
%    c_approx = c_approx + dt*duh2;

%% Hacemos el metodo para Diferencias de segundo orden (Método explicito)

  %  for i = 2:m-1
  %      for j = 2:n-1
  %          c_approx(i,j) = c_approx(i,j) + dt*((v/dx^2)*(c_approx(i+1,j) - 2*c_approx(i,j) + c_approx(i-1,j)) + (v/dy^2)*(c_approx(i,j+1) - 2*c_approx(i,j) + c_approx(i,j-1))); 
  %      end
  %  end


%% Imprimir solucion
    if k==0.01   %t/2
        %figure
        u_approx = reshape(c_approx,n,m);
        disp(u_approx)
        disp(c_approx)
        subplot(1,3,2)
        surf(x,y,u_approx)
        xlabel('x')
        ylabel('y')
        zlabel('c_approx')
    elseif k == 0.02   %t
        %figure
        u_approx = reshape(c_approx,n,m);
        %position=((k*10)/2)+1;
        disp(u_approx)
        disp(c_approx)
        subplot(1,3,3)
        surf(x,y,u_approx)
        xlabel('x')
        ylabel('y')
        zlabel('c_approx')
    end
    %disp(k)
    

end
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


end
