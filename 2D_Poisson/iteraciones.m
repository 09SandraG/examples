function [iterac,tx,tm,sdt] = iteraciones(n)
%% Evaluacion del tiempo para la funcion de poisson en 2D 
tx = zeros(n,1);
iterac = zeros(n,1);

tiemp = zeros(n,20);  %% Para determinar el tiempo de ejecucion
sd = 1;

while sd < 21
    for i = 1:n
        m = 10*i;
        [phi_approx,phi_exacta,x,y,tiempo,iter] = Poisson2D2(m,m,@phi,@f);
        iterac(i) = iter;
        tx(i) = m*m;
        tiemp(i,sd) = tiempo;
    end    
    sd = sd+1;
end
tm = mean(tiemp,2);
%disp(tm)
sdt = std(tiemp,0,2);
%disp(sdt)

%% Grafica del num de iteraciones vs tamaño del problema
%subplot(1,2,1)
figure
plot(tx,iterac);
title('Iteraciones del método de Gauss-Seidel para la ecuación de Poisson 2D');
xlabel('tamaño del problema');
ylabel('No. de iteraciones');

%% Grafica del tiempo de ejecución
%subplot(1,2,2)
figure
plot(tx,tm);
title('Tiempo de ejecución de la ecuación de Poisson 2D');
xlabel('tamaño del problema');
ylabel('segundos');
end