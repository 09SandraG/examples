function [g] = g(x,y,f)
% funcion que representa la invasion del tumor
% y el tratamiento en caso de tener alguno.
% Esta funcion depende de la densidad del tumor, por lo cual
% se manda a llamar a f.
% Aquí es donde se cambian o elegen la funcion de crecimiento y de
% tratamiento
%% Sin crecimiento
gc = 0;

%% crecimiento exponencial
% r>0: tasa de crecimiento instantanea
%r = 0.00012;
%gc = r*f;

%% Crecimiento logístico
% r>0: tasa de crecimiento instantanea
% K: capacidad de carga
%r = 0.0002;
%K = 10000;
%gc = r*f*(1-(f/K));

%% Crecimiento de Gompertz
% a: constante relacionada con la capacidad proliferativa de las celulas
% K: capacidad de carga
%a = 0.0003;
%K = 5000;
%gc = a*log(K/f)*f;

%% Tipos de tratamientos

%% Definicion de la función g
g = gc;
end