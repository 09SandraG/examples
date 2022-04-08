function [t,tx] = MedTiempo(n)
%% Evaluacion del tiempo para la funcion de poisson en 2D 
tx = zeros(n,1);
t = zeros(n,30);
sd = 1;

while sd < 31
    for i = 1:n
        m = 10*i;
        [phi_approx,phi_exacta,x,y,tiempo,cont11] = Poisson2D2(m,m,@phi,@f);
        t(i,sd) = tiempo;
        tx(i) = m;
    end    
    sd = sd+1;
end
tm = mean(t,2);
disp(tm)
sdt = std(t,0,2);
disp(sdt)



%% Grafica del tiempo de ejecución
%subplot(1,3,3)
plot(tx,tm);
title('Tiempo de ejecución de la ecuación de Poisson 2D');
xlabel('tamaño de la malla por lado');
ylabel('segundos');
end