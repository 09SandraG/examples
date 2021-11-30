function [t,tx] = MedTiempo(n)
%% Evaluacion del tiempo para la funcion de poisson en 3D

tx = zeros(n,1);
t = zeros(n,1);
for i = 1:n
    m = 10*i;
    [c_approx,x,y,tiempo] = ProliferacionInvasion2D(m,m,10,@f,@g,0.2);
    t(i) = tiempo;
    tx(i) = m;
end

%% Grafica del tiempo de ejecucion
plot(tx,t);
title('Tiempo de ejecucion del modelo Proliferación-invasión 2D');
xlabel('tamaño de la malla por lado');
ylabel('segundos');
end