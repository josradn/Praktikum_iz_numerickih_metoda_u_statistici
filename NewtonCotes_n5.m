%Zadatak 27 - Josipa Radnić, 1191240361
%Newton-Cotes n=5
function [I]=NewtonCotes_n5(a,b)
  n=5;
  h=(b-a)/n;
  x=linspace(a,b,n+1);
  I=5*h/288*(19*f(x(1))+75*f(x(2))+50*f(x(3))+50*f(x(4))+75*f(x(5))+19*f(x(6))); %netocno u literaturi
  I=2/sqrt(pi)*I;
end
