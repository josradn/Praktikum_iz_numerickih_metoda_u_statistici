%Zadatak 27 - Josipa Radnić, 1191240361
%Adaptivna Gaussova integracija
%u xi su pozitivne nultocke Legendreovih polinoma (integracijske tocke za N=2,4,6,8,10,12 iz tablice)
function [I,N]=Gaussova_adap(a,b,tol)
   h=(b-a)/2;
   m=(a+b)/2;
   n=1;

   Ans=Gauss(1,h,m);
   do
     Lastans=Ans;
     n=n+1;
     Ans=Gauss(n,h,m);
   until (abs(Ans-Lastans)<tol || n>6)

   if (n>6)
     flag=1
     return
   end

   I=2/sqrt(pi)*Ans;
   N=2*n;
end
