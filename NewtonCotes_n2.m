%Zadatak 27 - Josipa Radnić, 1191240361
%Simpsonova formula, n=2
function [I]=NewtonCotes_n2(a,b)
  I=2/sqrt(pi)*(b-a)/6*(f(a)+4*f((a+b)/2)+f(b));
end
