function [t,tx,tm,sdt] = iteraciones(n)
%% Evaluacion del tiempo y número de iteraciones aplicando el método de 
%Gauss-Seidel para la funcion de proliferación invasión
% n es el numero máximo de divisiones de lado de la malla
tx = zeros(n,1);
t = zeros(n,1);
tiemp = zeros(n,1);  %% Para determinar el tiempo de ejecucion
%sd = 1;

%while sd < 21  %%Numero de iteraciones para tomar promedio de los tiempos
    for i = 1:n
        m = 10*i;
        [tiempo,iter] = ProliferacionInvasion2D_MetodoImplicito(m,m);
        t(i) = iter;
        tx(i) = m;
        tiemp(i) = tiempo;
    end    
%    sd = sd+1;
%end
%tm = mean(tiemp,2);  %% Promedio de los tiempos de ejecución
tm = tiemp;
%disp(tm)
%sdt = std(tiemp,0,2);
sdt=1;
%disp(sdt)

%% Grafica del num de iteraciones vs tamaño del problema
%subplot(1,2,1)
%figure
%plot(tx,t);
%title('Iteraciones del método de Gauss-Seidel para el modelo de proliferación-Invasión');
%xlabel('tamaño de la malla por lado');
%ylabel('No. de iteraciones');

%% Grafica del tiempo de ejecución
%subplot(1,2,2)
figure
plot(tx,tm);
title('Tiempo de ejecución de la ecuación de Poisson 2D');
xlabel('tamaño de la malla por lado');
ylabel('segundos');
end