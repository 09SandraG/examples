function [xk] = gaussSeidel(A,b,x0,n,tol)

%variables de entrada

% A   :matriz de coeficientes del sistema A*x=b
% b   :vector de terminos independientes
%x0   :vector inicial
%tol  :toloerancia

%otras variables 

%n   :numero de variables
%ER  :norma de la diferencia de dos iteraciones consecutivas 
%xk  :k-esima iteracion

%fprintf('metodo de Gauss- Seidel:\n\n');
s='x';
%fprintf(' i');
%for i=1:n
%    fprintf('%12.1s %.0f',s,i);
%end
%fprintf('      ER\n');
ER=2; con =0;
while ER>tol
    con=con+1;
    for i=1:n
        ja=0; jd = 0;
        for j=1:i-1
            ja= ja +A(i,j)*xk(j);
        end
        for j=i+1:n
            jd= jd +A(i,j)*x0(j);
        end
        xk(i)=(b(i)-ja-jd)/A(i,i);
    end
    sum =0;
    for i=1:n
        sum= sum +(xk(i)-x0(i))*(xk(i)-x0(i));
    end
    ER= sqrt(sum);
 %   fprintf('%3.0f',con);
 %   for i=1:n
 %       fprintf('%15.6f',xk(i));
 %   end
 %  fprintf('%12.3e\n',ER);
    x0=xk;
    if con > 100
        fprintf('se excedio limite de iteraciones');
        ER=0;
    end
end
fprintf('Error: %12.3e\n',ER);
fprintf('Número de iteraciones: %3.0f \n',con);
end