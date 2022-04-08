function [it_total,tx,tiempo_Prom,sdt,error_Prom] = tiempos(n)
%% Evaluacion del tiempo y número de iteraciones aplicando el método de 
%Gauss-Seidel para la funcion de proliferación invasión
% n es el numero máximo de divisiones de lado de la malla
tx = zeros(n,1);  %Tamaño del problema o malla
tiemp = zeros(5,1);  %% Para determinar el tiempo de ejecucion
aux_iter = zeros(5,1);
aux_error = zeros(5,1);
tiempo_Prom = zeros(n,1);
it_total = zeros(n,1);   %Número de iteraciones (para métodos iterativos)
error_Prom = zeros(n,1);
sdt = zeros(n,1);

for i = 1:n
    sd = 1;
    m = 10*i;
    disp(m*m)
    while sd < 6  %%Numero de iteraciones para tomar promedio de los tiempos
        disp(sd)
        [tiempo,iter,error] = ProliferacionInvasion2D_MetodoImplicito(m,m);
        tiemp(sd) = tiempo;
        aux_iter(sd) = iter;
        aux_error(sd) = error;
        sd = sd+1;
    end
    it_total(i) = mode(aux_iter); 
    tx(i) = m*m;
    tiempo_Prom(i) = mean(tiemp); %% Guardamos tiempos Promedio de los tiempos de ejecución
    error_Prom(i) = mean(aux_error);
    sdt(i) = std(tiemp);
end     
%tm = tiemp;
%disp(tm)
%sdt=1;
%disp(sdt)

%% Grafica de tamaño del problema vs numero de iteraciones
%subplot(1,2,1)
figure
plot(tx,it_total);
axis([100 inf 2 7]);
title('Iterations. Gauss-Seidel method. Proliferation-invasion Model 2D');
xlabel('problem size');
ylabel('number of iterations');

%% Grafica del tamaño del problema vs tiempo de ejecución 
%subplot(1,2,2)
figure
plot(tx,tiempo_Prom);
title('Solution times. Proliferation-invasion Model 2D');
xlabel('problem size');
ylabel('seconds');

%% Grafica del num de iteraciones vs tolerancia o error
%subplot(1,2,1)
%figure
%plot(it_total,error_Prom);
%title('Iterations vs tolerance. Gauss-Seidel method. Proliferation-invasion Model 2D');
%xlabel('problem size');
%ylabel('number of iterations');

end