function [t,tx] = MedTiempo(n)
%% Evaluacion del tiempo para la funcion de poisson en 3D
tx = zeros(n,1);
t = zeros(n,20);
sd = 1;

while sd < 21
    for i = 1:n
        m = 10*i;
        [phi_approx, phi_exacta,x,y,z,tiempo,cont11] = Poisson3D2(m,m,m,@phi,@f);
        t(i,sd) = tiempo;
        tx(i) = m;
    end
    sd = sd+1;
end

tm = mean(t,2);
disp(tm)
sdt = std(t,0,2);
disp(sdt)


%% Grafica del tiempo de ejecucion
plot(tx,tm);
title('Tiempo de ejecucion de la ecuación de Poisson 3D');
xlabel('tamaño de la malla por lado');
ylabel('segundos');
end