%Zadatak 27 - Josipa Radnić, 1191240361
%Milenovo ili Rombergovo pravilo, n=4
function [I]=NewtonCotes_n4(a,b)
  n=4;
  h=(b-a)/n;
  x=linspace(a,b,n+1);
  I=2*h/45*(7*f(x(1))+32*f(x(2))+12*f(x(3))+32*f(x(4))+7*f(x(5)));
  I=2/sqrt(pi)*I;
end
