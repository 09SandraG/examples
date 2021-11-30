function [t,tx] = MedTiempo(n)
%% Evaluacion del tiempo para la funcion de poisson en 2D 
tx = zeros(n,1);
t = zeros(n,1);
for i = 1:n
    m = 10*i;
    [phi_approx,phi_exacta,x,y,tiempo] = Poisson2D2(m,m,@phi,@f);
    t(i) = tiempo;
    tx(i) = m;
end

%% Grafica del tiempo de ejecución
%subplot(1,3,3)
plot(tx,t);
title('Tiempo de ejecución de la ecuación de Poisson 2D');
xlabel('tamaño de la malla por lado');
ylabel('segundos');
end