function [t,tx] = MedTiempo(n)
%% Evaluacion del tiempo para la funcion de poisson en 3D

tx = zeros(n,1);
t = zeros(n,1);
for i = 1:n
    m = 10*i;
    t1t = @() Poisson3D2(m,m,m,@phi,@f);
    t(i) = timeit(t1t);
    tx(i) = m;
end

%% Grafica del tiempo de ejecucion
plot(tx,t);
end