%Zadatak 27 - Josipa Radnić, 1191240361
%kompozitna trapezna formula, n=1
function [I]=kompozitna_trapezna(a,b,N)
  h=(b-a)/N;
  sum=(f(a)+f(b))/2;
  for i=1:(N-1)
    sum=sum+f(a+i*h);
  end
  I=2/sqrt(pi)*h*sum;
end


