%Zadatak 27 - Josipa Radnić, 1191240361
%Newton-Cotes n=8
function [I]=NewtonCotes_n8(a,b)
  n=8;
  h=(b-a)/n;
  x=linspace(a,b,n+1);
  I=4*h/14175*(989*f(x(1))+5888*f(x(2))-928*f(x(3))+10496*f(x(4))-4540*f(x(5))+10496*f(x(6))-928*f(x(7))+5888*f(x(8))+989*f(x(9)));
  I=2/sqrt(pi)*I;
end
