%Zadatak 27 - Josipa Radnić, 1191240361
%Newtonovo 3/8 pravilo, n=3
function [I]=NewtonCotes_n3(a,b)
  n=3;
  h=(b-a)/n;
  x=linspace(a,b,n+1);
  I=3*h/8*(f(x(1))+3*f(x(2))+3*f(x(3))+f(x(4)));
  I=2/sqrt(pi)*I;
end

