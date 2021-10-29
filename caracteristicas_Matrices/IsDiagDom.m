function [isdom] = IsDiagDom( A ) 
isdom = true;
for r = 1:size(A,1)
    rowdom = 2 * abs(A(r,r)) > sum(abs(A(r,:)));
    isdom = isdom && rowdom;
end
%if isdom == 0
%    disp (['Matrix A is not diagonally-dominant']);
%elseif isdom == 1
%        disp (['Matrix A is diagonally-dominant']); 
%end
end