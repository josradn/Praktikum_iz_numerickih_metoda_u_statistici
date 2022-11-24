%Zadatak 27 - Josipa Radnić, 1191240361
%kompozitna Simpsonova formula, n=2
function [I]=kompozitna_Simpsonova(a,b,N)
  h=(b-a)/(2*N);
  sum=f(a)+f(b);
  for i=1:(N-1)
    sum=sum+4*f(a+(2*i-1)*h)+2*f(a+(2*i)*h);
  end
  I=2/sqrt(pi)*h/3*sum;
end


