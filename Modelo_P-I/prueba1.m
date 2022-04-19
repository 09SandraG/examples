
A = [1 2 3 3 4 1 6 4 9 1; 8 0 2 3 0 1 9 4 1 1; 1 1 3 2 4 3 2 4 5 9; 0 2 4 3 9 1 4 4 2 0; 1 2 3 2 4 8 6 4 9 8; 8 0 3 5 4 0 6 4 0 0; 3 3 3 3 6 1 7 9 9 1; 5 2 1 2 3 6 6 4 5 9; 0 2 5 4 4 9 2 4 5 1; 1 2 5 0 4 1 8 4 4 3];
B = [3 6 0 2 7 1 9 6 0 6; 3 4 6 2 7 1 9 6 0 5; 2 6 8 9 3 1 6 6 0 9; 1 2 1 3 7 6 8 6 0 7; 5 2 1 8 7 7 2 6 8 5; 2 6 8 2 8 9 0 4 0 3; 3 5 0 3 7 9 7 5 0 2; 2 4 5 2 6 1 2 6 0 0; 3 4 1 3 8 1 9 6 0 8; 0 2 5 4 4 9 2 4 5 1];
    
[C1, tiempo1] = mult1(A,B);

tStart2 = cputime;
C2 = A*B;
tiempo2 = cputime - tStart2;

disp(C1)
disp(tiempo1)
disp(C2)
disp(tiempo2)

function [m,n,q] = check(A,B)
    [m,n] = size(A);
    [p,q] = size(B);
    assert(n == p, 'Las matrices no son compatibles');
end

function [C,tiempo] = mult1(A,B)
    [m,n,q] = check(A,B);
    C = zeros(m,q);
    tStart = cputime;
    for i=1:m
        for j=1:q
            for k=1:n
                C(i,j) = C(i,j) + A(i,k)*B(k,j);
            end
        end
    end
    tiempo = cputime - tStart;
end

