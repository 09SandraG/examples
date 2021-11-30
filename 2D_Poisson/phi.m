function [phi] = phi(x,y)
%phi = 2.*exp(2*x+y);
alp = 1.0;
phii = 0.0;
phi=tanh(1-alp*(x*tan(phii)-y));
end