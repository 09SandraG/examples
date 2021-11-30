function [f] = f(x,y)
%f = 10.*exp(2*x+y); 
alp = 1.0;
phii = 0.0;
f = 2*alp^2*tanh(alp*(y - x*tan(phii)) + 1)*tan(phii)^2*(tanh(alp*(y - x*tan(phii)) + 1)^2 - 1) + 2*alp^2*tanh(alp*(y - x*tan(phii)) + 1)*(tanh(alp*(y - x*tan(phii)) + 1)^2 - 1);
end