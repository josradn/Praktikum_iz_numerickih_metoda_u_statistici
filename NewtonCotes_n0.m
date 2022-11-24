%Zadatak 27 - Josipa Radnić, 1191240361
%pravokutna formula, n=0
function [I]=NewtonCotes_n0(a,b,x)
  I=2/sqrt(pi)*f(x)*(b-a);
end

%probat cemo za razlicite x-eve, najcesci su x=a, x=b i x=(a+b)/2
