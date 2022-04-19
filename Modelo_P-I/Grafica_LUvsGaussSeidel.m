%%%% Grafica de los tiempos de ejecuciones del método de  
%descomposición LU vs Gauss-Seidel

x = [100 400 900 1600 2500 3600 4900 6400 8100 10000 12100 14400 16900];   % Tamaños del problema

y1 = [0 0.0008 0.0088 0.0387 0.1260 0.3214 0.6435 1.1753 2.2112 4.5456 8.3655 18.4929 32.5912];  % Datos de los tiempos de descomposición LU
y2 = [0.0107 0.0229 0.0364 0.0551 0.0883 0.1214 0.1743 0.2596 0.4191 0.8221 1.2545 2.5424 3.2981];  % Datos de los tiempos de Gauss-Seidel

figure
plot(x,y1,'r-.',x,y2,'b--','LineWidth',2);
legend({'times LU descomposition','times Gauss-Seidel'},'location','northwest');
