%Zadatak 27 - Josipa Radnić, 1191240361
%Rombergov algoritam
function [I,N]=Rombergov_alg(a,b,tol,MaxLevel)
  LastRow=zeros(1,MaxLevel+1);
  CurrRow=zeros(1,MaxLevel+1);
  Level=0;
  N=1;
  h=(b-a)/N;
  %CurrRow(0)=T(f,a,b,1)
  CurrRow(1)=1/2*h*(f(a)+f(b));

  do
    N=2*N;
    LastRow=CurrRow;
    sum=0;
    h=(b-a)/N;
    for i=1:(N/2)
      sum=sum+f(a+(2*i-1)*h);
    end
    CurrRow(1)=1/2*LastRow(1)+h*sum;
    Level=Level+1;
    Coeff=1;
    for j=1:Level
      Coeff=4*Coeff;
      CurrRow(j+1)=(Coeff*CurrRow(j)-LastRow(j))/(Coeff-1);
    end
  until (abs(CurrRow(Level+1)-LastRow(Level))<tol || Level>MaxLevel)

  if(Level>MaxLevel)
    flag=1
    return
  end

  I=2/sqrt(pi)*CurrRow(Level+1);
end
