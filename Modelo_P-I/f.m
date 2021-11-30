function [f] = f(x,y,t)
% f es la funcion de densidad que depende de x,y o podemos dar una constante. 
% Esta función se llama en g pues el crecimiento depende de la densidad.

f = 0.01 + (x*y) + 0*t;
end