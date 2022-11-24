%Zadatak 27 - Josipa Radnić, 1191240361
format long e

a=0;
b=2;
x_left=a;
x_right=b;
x_middle=(a+b)/2;
n=200;
tol1=5e-4; %u materijalima koristena tolerancija
tol2=1e-6;

%Usporediti cemo pravu vrijednost funkcije za x=2 s vrijednostima dobivenima Newton-Cotesovim metodama
y_real=erf(2)
[I0_left]=NewtonCotes_n0(a,b,x_left)
[I0_right]=NewtonCotes_n0(a,b,x_right)
[I0_middle]=NewtonCotes_n0(a,b,x_middle)
[I1]=NewtonCotes_n1(a,b)
[I2]=NewtonCotes_n2(a,b)
[I3]=NewtonCotes_n3(a,b)
[I4]=NewtonCotes_n4(a,b)
[I5]=NewtonCotes_n5(a,b)
[I8]=NewtonCotes_n8(a,b)
%y_real = 9.953222650189527e-01
%I0_left = 2.256758334191025e+00
%I0_right = 4.133397070818411e-02
%I0_middle = 8.302149948411894e-01
%I1 = 1.149046152449605e+00
%I2 = 9.364920473773278e-01
%I3 = 9.729158106527891e-01
%I4 = 9.989205514952164e-01
%I5 = 9.972650025700477e-01
%I8 = 9.953279462514626e-01
%Vidimo da kako povecavamo n da se vrijednost integrala priblizava pravoj vrijednosti integrala
%Takoder vidimo da za pravokutnu formulu srednja vrijednost najbolje daje aproksimaciju od svih 3


%Prikazimo graficki funkciju na intervalu [0,2]:
figure(1)
x=linspace(a,b,n);
plot(x,erf(x),'red');

%Sada cemo zajedno s funkcijom graficki prikazati Newton-Cotesove metode:
for i=1:n
  y0_left(i)=NewtonCotes_n0(a,x(i),a);
  y0_right(i)=NewtonCotes_n0(a,x(i),x(i));
  y0_middle(i)=NewtonCotes_n0(a,x(i),(x(i)+a)/2);
  y1(i)=NewtonCotes_n1(a,x(i));
  y2(i)=NewtonCotes_n2(a,x(i));
  y3(i)=NewtonCotes_n3(a,x(i));
  y4(i)=NewtonCotes_n4(a,x(i));
  y5(i)=NewtonCotes_n5(a,x(i));
  y8(i)=NewtonCotes_n8(a,x(i));
end

figure(2)
plot(x,erf(x),'red');
hold on;
plot(x,y0_left,'green');
hold on;
plot(x,y0_right,'blue');
hold on;
plot(x,y0_middle,'cyan');
title("Grafovi za Newton Cotesove pravokutne formule (n=0) uz prikaz funkcije");
h=legend("funkcija","left", "right", "middle");
legend(h,"location","northeastoutside");
hold off;

figure(3)
plot(x,erf(x),'red');
hold on;
plot(x,y1,'green');
hold on;
plot(x,y2,'blue');
hold on;
plot(x,y3,'cyan');
hold on;
plot(x,y4,'magenta');
hold on;
plot(x,y5,'yellow');
hold on;
plot(x,y8,'black');
title("Grafovi za Newton Cotesove formule uz prikaz funkcije");
h=legend("funkcija","n=1 (trapezna)", "n=2 (Simpsonova)", "n=3 (Newtonovo 3/8 pravilo)", "n=4 (Rombergovo pravilo)", "n=5","n=8");
legend(h,"location","northeastoutside");
hold off;


%Sada cemo usporediti kompozitne Newton-Cotesove formule s funkcijom
%Prvo pronadimo koji bi broj N trebao biti preko ocjene greske za odredenu toleranciju:

for i=1:n
  d2fvrj(i)=abs(d2f(x(i)));
  d4fvrj(i)=abs(d4f(x(i)));
end
maxd2f=max(d2fvrj);
maxd4f=max(d4fvrj);

%TOL1
N_trap=sqrt((b-a)^3/(tol1*12)*maxd2f);
N_Sim=nthroot(((b-a)^5/(tol1*180)*maxd4f),4);
%N_trap = 5.163977794943222e+01 -> N_trap>=52
%N_Sim = 8.082062018706493e+00 -> N_Sim>=9
for i=2:52
  I_komp_trap_tol1(i)=kompozitna_trapezna(a,b,i);
end
I_komp_trap_tol1;
I_komp_trap_tol1_zadnja=I_komp_trap_tol1(52)
%Prava vrijednost od integrala je 9.953222650189527e-01, dok zadnja iteracija
%u kompozitnoj trapeznoj formuli daje 9.953120767103369e-01, sto je blizu pravoj vrijednosti
for i=2:9
  I_komp_Sim_tol1(i)=kompozitna_Simpsonova(a,b,i);
end
I_komp_Sim_tol1;
I_komp_Sim_tol1_zadnja=I_komp_Sim_tol1(9)
%Prava vrijednost od integrala je 9.953222650189527e-01, dok zadnja iteracija
%u kompozitnoj Simpsonovoj formuli daje  9.906049384180480e-01, sto bas nije idealno

%TOL2
N_trap=sqrt((b-a)^3/(tol2*12)*maxd2f);
N_Sim=nthroot(((b-a)^5/(tol2*180)*maxd4f),4);
%N_trap = 1.154700538379251e+03 -> N_trap>=1155
%N_Sim = 3.821771168817467e+01 -> N_Sim>=39
for i=2:1155
  I_komp_trap_tol2(i)=kompozitna_trapezna(a,b,i);
end
I_komp_trap_tol2;
I_komp_trap_tol2_zadnja=I_komp_trap_tol2(1155)
%Prava vrijednost od integrala je 9.953222650189527e-01, dok zadnja iteracija
%u kompozitnoj trapeznoj formuli daje 9.953222443626926e-01, sto je jako blizu pravoj vrijednosti
for i=2:39
  I_komp_Sim_tol2(i)=kompozitna_Simpsonova(a,b,i);
end
I_komp_Sim_tol2;
I_komp_Sim_tol2_zadnja=I_komp_Sim_tol2(39)
%Prava vrijednost od integrala je 9.953222650189527e-01, dok zadnja iteracija
%u kompozitnoj Simpsonovoj formuli daje  9.945398992238774e-01, sto je blizu pravoj vrijednosti

%Mozemo primjetiti da kompozitna trapezna formula daje bolju tocnost integrala, dok je potrebno
%vise podintervala napraviti kako bi se to postiglo. No, kompozitna Simpsonova za tol1 daje
%manju tocnost, no za puno manji broj podintervala. Tako da zapravo bi ocjena greske kod
%kompozitne trapezne formule bila manja, no broj potrebnih podintervala veci.
%U slucaju tol2 kompozitna Simpsonova trazi puno manji broj podintervala, ali je svejedno
%dao zadovoljavajuci rezultat integrala. Dakle, ovisi i o toleranciji, sto je tolerancija manja
%(dakle ostrija) to ce kompozitna Simpsonova dat bolji rezultat za puno manji broj podintervala
%nego kompozitna trapezna.


%Sada cemo za obe tolerancije pokazati graficki kompozitne formule s funkcijom
for i=1:n
  d2fvrj=0;
  d4fvrj=0;
  maxd2f=0;
  maxd4f=0;
  for j=1:i
    d2fvrj(j)=abs(d2f(x(j)));
    d4fvrj(j)=abs(d4f(x(j)));
  end
  maxd2f=max(d2fvrj);
  maxd4f=max(d4fvrj);
  N_trap=0;
  N_Sim=0;
  N_trap=sqrt((x(i)-a)^3/(tol1*12)*maxd2f);
  N_Sim=nthroot(((x(i)-a)^5/(tol1*180)*maxd4f),4);
  y_komp_trap_tol1(i)=kompozitna_trapezna(a,x(i),ceil(N_trap));
  y_komp_Sim_tol1(i)=kompozitna_Simpsonova(a,x(i),ceil(N_Sim));
end

for i=1:n
  d2fvrj=0;
  d4fvrj=0;
  maxd2f=0;
  maxd4f=0;
  for j=1:i
    d2fvrj(j)=abs(d2f(x(j)));
    d4fvrj(j)=abs(d4f(x(j)));
  end
  maxd2f=max(d2fvrj);
  maxd4f=max(d4fvrj);
  N_trap=0;
  N_Sim=0;
  N_trap=sqrt((x(i)-a)^3/(tol2*12)*maxd2f);
  N_Sim=nthroot(((x(i)-a)^5/(tol2*180)*maxd4f),4);
  y_komp_trap_tol2(i)=kompozitna_trapezna(a,x(i),ceil(N_trap));
  y_komp_Sim_tol2(i)=kompozitna_Simpsonova(a,x(i),ceil(N_Sim));
end

figure(4)
plot(x,erf(x),'red');
hold on;
plot(x,y_komp_trap_tol1, 'green');
hold on;
plot(x,y_komp_Sim_tol1, 'blue');
title("Grafovi za kompozitne Newton Cotesove formule uz prikaz funkcije (tol=5e-4)");
h=legend("funkcija", "kompozitna trapezna","kompozitna Simpsonova");
legend(h,"location","northeastoutside");
hold off;

figure(5)
plot(x,erf(x),'red');
hold on;
plot(x,y_komp_trap_tol2, 'green');
hold on;
plot(x,y_komp_Sim_tol2, 'blue');
title("Grafovi za kompozitne Newton Cotesove formule uz prikaz funkcije (tol=1e-6)");
h=legend("funkcija", "kompozitna trapezna","kompozitna Simpsonova");
legend(h,"location","northeastoutside");
hold off;
%Vidimo da kod tol=1e-6 kompozitna Simpsonova je bliza nego kod tol=5e-4, no i dalje je
%veca razlika nego kod trapezne. Za obe tolerancije kompozitna trapezna preklapa s funkcijom
%te je ona dosta bolja sto se tice tocnosti.


%Sada cemo usporediti Newton-Cotesove formule za n=5 i n=8 i kompozitne Newton-Cotesove formule za tol=1e-6:
figure(6)
plot(x,erf(x),'red');
hold on;
plot(x,y5,'yellow');
hold on;
plot(x,y8,'black');
hold on;
plot(x,y_komp_trap_tol2, 'green');
hold on;
plot(x,y_komp_Sim_tol2, 'blue');
title("Grafovi za Newton-Cotesove i kompozitne Newton-Cotesove formule uz prikaz funkcije (tol=1e-6)");
h=legend("funkcija", "n=5", "n=8", "kompozitna trapezna","kompozitna Simpsonova");
legend(h,"location","northeastoutside");
hold off;

%Vidimo da za n=5,n=8 i kompozitnu trapeznu jako poklapaju s funkcijom i kompozitna Simpsonova je blizu njima


%Sada pogledajmo Gaussovu integraciju i adaptivnu Gaussovu integraciju

for i=1:5
  [I_Gauss(i),~]=Gaussova_int(a,b,i);
end
%I_Gauss=8.085190034231082e-01, 7.379688843032787e-01, 1.052847181681450e+00,
%        9.594040243551228e-01, 1.084526582069845e+00
%Rezultati nisu najidealniji
greske_Gauss=abs(I_Gauss-y_real)
%greske_Gauss=1.868032615958445e-01, 2.573533807156740e-01, 5.752491666249726e-02,
%             3.591824066382987e-02, 8.920431705089260e-02
%Vidimo da za i=4(N=5) najmanju gresku napravi

tic
[I_Gauss_adap_tol1,N_Gauss_adap_tol1]=Gaussova_adap(a,b,tol1)
toc
%I_Gauss_adap_tol1 = 9.953222250622672e-01
%N_Gauss_adap_tol1 = 6.000000000000000e+00
tic
[I_Gauss_adap_tol2,N_Gauss_adap_tol2]=Gaussova_adap(a,b,tol2)
toc
%I_Gauss_adap_tol2 = 9.953222646324572e-01
%N_Gauss_adap_tol2 = 8.000000000000000e+00
greske_Gauss_adap_tol1=abs(I_Gauss_adap_tol1-y_real)
greske_Gauss_adap_tol2=abs(I_Gauss_adap_tol2-y_real)
%greske_Gauss_adap_tol1 = 3.995668551581844e-08
%greske_Gauss_adap_tol2 = 3.864955022692129e-10
%Greske su ocekivano puno manje u adaptivnoj metodi


%Sada cemo graficki usporediti Gaussove integracijske formule s funkcijom
for i=1:n
  [y_Gauss(i),~]=Gaussova_int(a,x(i),4);
  [y_Gauss_adap_tol1(i),~]=Gaussova_adap(a,x(i),tol1);
  [y_Gauss_adap_tol2(i),~]=Gaussova_adap(a,x(i),tol2);
end

figure(7)
plot(x,erf(x),'red');
hold on;
plot(x,y_Gauss_adap_tol1,'green');
hold on;
plot(x,y_Gauss_adap_tol2,'blue');
hold on;
plot(x,y_Gauss, 'yellow');
title("Grafovi za adaptivnu Gaussovu integraciju i Gaussovu uz prikaz funkcije");
h=legend("funkcija", "adaptivna Gaussova tol=5e-4","adaptivna Gaussova tol=1e-6","Gauss");
legend(h,"location","northeastoutside");
hold off;
%Vidimo i na grafu da se adaptivna Gaussova za obe tolerancije jako dobro preklapa s funkcijom
%dok obicna Gaussova odskace od funkcije. Koristiti cemo adaptivnu Gaussovu integraciju za daljnje usporedbe


%Jos cemo pogledati Rombergov algoritam i usporediti ga s adaptivnim Gaussovim
%Za tol1 MaxLevel manji od 4 vraca gresku, te za tol2 MaxLevel manji od 5 vraca gresku, tako su izabrani brojevi
tic
[I_Romb_tol1,N_Romb_tol1]=Rombergov_alg(a,b,tol1,4)
toc
%I_Romb_tol1 = 9.953224646647567e-01
%N_Romb_tol1 = 1.600000000000000e+01
tic
[I_Romb_tol2,N_Romb_tol2]=Rombergov_alg(a,b,tol2,5)
toc
%I_Romb_tol2 = 9.953222649586621e-01
%N_Romb_tol2 = 3.200000000000000e+01
greske_Romb_tol1=abs(I_Romb_tol1-y_real)
greske_Romb_tol2=abs(I_Romb_tol2-y_real)
%greske_Romb_tol1 = 1.996458039998927e-07
%greske_Romb_tol2 = 6.029066135226913e-11
%Vidimo da je jako blizu pravoj vrijednosti integrala i da su greske jako male

%Sada cemo graficki usporediti Rombergov algoritam i adaptivnu Gaussovu integraciju s funkcijom za obe tolerancije
for i=1:n
  [y_Romb_tol1(i),~]=Rombergov_alg(a,x(i),tol1,4);
  [y_Romb_tol2(i),~]=Rombergov_alg(a,x(i),tol2,5);
end

figure(8)
plot(x,erf(x),'red');
hold on;
plot(x,y_Gauss_adap_tol1,'green');
hold on;
plot(x,y_Gauss_adap_tol2,'blue');
hold on;
plot(x,y_Romb_tol1,'cyan');
hold on;
plot(x,y_Romb_tol2,'magenta');
title("Grafovi za adaptivnu Gaussovu integraciju i Rombergov algoritam uz prikaz funkcije");
h=legend("funkcija", "adaptivna Gaussova tol=5e-4","adaptivna Gaussova tol=1e-6", "Rombergov tol=5e-4", "Rombergov tol=1e-6");
legend(h,"location","northeastoutside");
hold off;
%Vidimo da se linije za svaku toleranciju jako dobro preklapaju s funkcijom, tako da rezultati su zadovoljavajuci

%Jos cemo samo prikazati sve navedene metode na jednom grafu radi usporedbe za jednu toleraciju (tol2)
%Dakle: Newton-Cotesove formule za n=5 i n=8,kompozitna trapezna i Simpsonova, adaptivna Gaussova integracijske i Rombergov algoritam


figure(9)
plot(x,erf(x),'red');
hold on;
plot(x,y5,'green');
hold on;
plot(x,y8,'blue');
hold on;
plot(x,y_komp_trap_tol2, 'cyan');
hold on;
plot(x,y_komp_Sim_tol2, 'magenta');
hold on;
plot(x,y_Gauss_adap_tol2,'yellow');
hold on;
plot(x,y_Romb_tol2,'black');
title("Grafovi za različite metode uz prikaz funkcije (tol=1e-6)");
h=legend("funkcija", "Newton-Cotes n=5","Newton-Cotes n=8","kompozitna trapezna","kompozitna Simpsonova","adaptivna Gaussova", "Rombergov alg");
legend(h,"location","northeastoutside");
hold off;
%Vidimo da je kompozitna Simpsonova najlosija od svih drugih metoda, te da su sve druge jako dobre
%Kada bi uzeli u obzir brzinu izvrsavanja algoritma te preciznost i ocjenu greske algoritma,
%kompozitna trapezna dosta vremena treba za izvrsavanje sto nije najbolje rjesenje.
%Najmanje vremena za izvrsavanje te najmanja greska se dogada u adaptivnoj Gaussovoj integraciji
%i Rombergovom algoritmu. Po vremenu izvrsavanja sto sam gore testirala s funkcijama tic i toc
%moze se primjetiti da je adaptivna Gaussova integracija brza u izvsravanju nego Rombergov.


%Znamo da je integral od 0 do +besk od funkcije e^(-x^2) jednak sqrt(pi)/2 (dokazano na raznim kolegijima)
%pa je tada fi(beskonacno)=2/sqrt(pi)*(sqrt(pi)/2)=1
%Mozemo to provjeriti i pomocu funkcije u matlabu:
pkg load symbolic
beskonacno=erf(sym(Inf))
%Symbolic pkg v2.9.0: Python communication link active, SymPy v1.5.
%beskonacno = (sym) 1

%fi(+besk) se moze naci pomocu Gaussove formule i nultocaka Hermitskog polinoma.
%Hermitski polinomi su ortogonalni s obzirom na tezinsku fju w(x), tj. za hermitske polinome
%H_m i H_n gdje je m,n {0,1,2,..} vrijedi integral od -besk do +besk od H_m(x)*H_n(x)w(x)
%je jednak 0 ako je m!=n, odnosno sqrt(pi)*2^n*n!delta(n,m) gdje je w(x)=e^(-x^2), a delta je
%kroneckerova delta koja je 1 kad je m=n, inace 0.
%Kako je navedena tezinska funkcija zapravo podintegralna funkcija u ovom zadatku,
%vrijednost fi(+besk) je zapravo norma Hermitskog polinoma reda 0 pomnozena tezinskom funkcijom.
%Za n=0 imamo sqrt(pi)*2^0*0!*delta(0,0)=sqrt(pi)*1*1*1=sqrt(pi)
%Podintegralna funkcija je parna funkcija, time je integral od 0 do +besk pomnozen s 2
%zapravo jednak integralu od -besk do +besk od podintegralne funkcije.
%Kako mi imamo da je integral od -besk do +besk od podintegralne fje (e^(-x^2)) jednak sqrt(pi)
%tada je integral od 0 do +besk od 2/sqrt(pi)*e^(-x^2)=1/sqrt(pi)*(2*e^(-x^2)) jednak
% 1/sqrt(pi)*sqrt(pi)=1
erf(6) %vec daje zeljenu vrijednost
