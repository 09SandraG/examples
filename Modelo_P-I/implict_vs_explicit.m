function [tx,tiempo_Prom_explicit,tiempo_Prom_implicit,sdte,sdti] = implict_vs_explicit(n)
%%%% Comparación entre los tiempos totales 
% n: numero de iteraciones y aumento de la malla

tx = zeros(n,1);  %Tamaño del problema o malla
tiemp_explicit = zeros(5,1);  %% Para determinar el tiempo de ejecucion explicito
tiemp_implicit = zeros(5,1);  %% Para determinar el tiempo de ejecucion implicito
tiempo_Prom_explicit = zeros(n,1);
tiempo_Prom_implicit = zeros(n,1);
sdte = zeros(n,1);
sdti = zeros(n,1);

for i = 1:n
    sd = 1;
    m = 10*i;
    disp(m*m)
     while sd < 6  %%Numero de iteraciones para tomar promedio de los tiempos
        disp(sd)
        [tiempoe] = ProliferacionInvasion2D_MetodoExplicito(m,m);
        tiemp_explicit(sd) = tiempoe;
        [tiempoi,iter,error] = ProliferacionInvasion2D_MetodoImplicito(m,m);
        tiemp_implicit(sd) = tiempoi;
        sd = sd+1;
     end
    tx(i) = m*m;
    tiempo_Prom_explicit(i) = mean(tiemp_explicit); %% Guardamos tiempos Promedio de los tiempos de ejecución
    tiempo_Prom_implicit(i) = mean(tiemp_implicit); %% Guardamos tiempos Promedio de los tiempos de ejecución
    sdte(i) = std(tiemp_explicit);
    sdti(i) = std(tiemp_implicit);
end     

figure
plot(tx,tiempo_Prom_explicit,'r-.',tx,tiempo_Prom_implicit,'b--','LineWidth',2);
legend({'times explicit method','times implicit method'},'location','northwest');
