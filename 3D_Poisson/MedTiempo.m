function [t,tx] = MedTiempo(n)
%% Evaluacion del tiempo para la funcion de poisson en 3D

tx = zeros(n,1);
t = zeros(n,1);
for i = 1:n
    m = 10*i;
    [phi_approx, phi_exacta,x,y,z,tiempo] = Poisson3D(m,m,m,@phi,@f);
    t(i) = tiempo;
    tx(i) = m;
end

%% Grafica del tiempo de ejecucion
plot(tx,t);
title('Tiempo de ejecucion');
xlabel('tamaño de la malla');
ylabel('segundos');
end