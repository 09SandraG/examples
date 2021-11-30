function [x_approx] = Jacobi(A,b,x0,n,tol,maxit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Métodod iterativo de Jacobi              %%
% Modificado por Sandra García. 28/10/2021 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Entradas %%%
% n        numero de ecuaciones e incognitas
% A        la matriz de nxn
% b        Vector de nx1
% x0       Vector inicial
% tol      tolerancia
% maxit    numero maximo de iteraciones

%%%%% Salidas %%%%%%
% x_approx Solucion aproximada

%%%% Inicio del programa %%
x_approx = zeros(n,1);
k = 1;
err = 10;
while (k<=maxit) && (err > tol)
    for i = 1:n
        if A(i,i) == 0
            err = 0;
            break;
        end
        sum1 = 0;
        for j = 1:n
            if i ~= j
                sum1 = sum1 + A(i,j)*x0(j);
            end
        end
        x_approx(i) = (-sum1 + b(i))/A(i,i);
    end
    %%% El error se calcula con la norma euclideana
    %fprintf('%3.0f',k);
    sum = 0;
    for i = 1:n
        sum = sum + (x_approx(i) - x0(i))*(x_approx(i) - x0(i));
    end
    err = sqrt(sum);
    %fprintf('%12.3e\n',err);
    k=k+1;
    x0 = x_approx;
end
fprintf('Despues de %3.0f iteraciones el error de la aproximación es: %3.6e\n',k-1,err);
%disp(k)
end