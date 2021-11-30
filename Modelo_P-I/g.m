function [g] = g(x,y,f,t)
% funcion que representa la invasion del tumor
% y el tratamiento en caso de tener alguno.
% Esta funcion depende de la densidad del tumor, por lo cual
% se manda a llamar a f.
% Aquí es donde se cambian o elegen la funcion de crecimiento y de
% tratamiento

%% crecimiento exponencial
% r>0: tasa de crecimiento instantanea
%r = 0.2;
%gc = r*f(x,y);

%% Crecimiento logístico
% r>0: tasa de crecimiento instantanea
% K: capacidad de carga
r = 0.02 + t*0;
K = 100;
gc = r*f(x,y,t)*(1-(f(x,y,t)/K));

%% Crecimiento de Gompertz
% a: constante relacionada con la capacidad proliferativa de las celulas
% K: capacidad de carga
%a = 0.3;
%K = 5;
%gc = a*log(K/f(x,y))*f(x,y);

%% Tipos de tratamientos

%% Definicion de la función g
g = gc;
end