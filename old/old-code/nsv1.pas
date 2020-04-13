program sigtcf;

uses lnrplr,crt, complex1;

label 77;
const pi=3.14159265;
      mev=0.511;
var d,df,aa,gv,z1,z2,a1,a2,mn,toch,m,tt:real;
    i,n,id:integer;
    an:char;
    send_ef:real;
    send_ei:real;
    calcul:integer;{для показа прогресса счета}
    f:text;

function kf(ef:real):real;
begin
kf:=sqrt(2*m*ef)
end;

function lf(ef:real):real;
begin
lf:=(z1+1)*z2*aa*sqrt(m/2/ef)
end;

function ki(ei:real):real;
begin
  ki:=sqrt(2*m*ei);
end;

function li(ei:real):real;
begin
  li:=z1*z2*aa*sqrt(m/2/ei);
end;

{function fun1(t:real):real;
var a,b,c,x,rez:complex;
    lli,llf,kki,kkf,qq,qqq,qqqq,per,ef,ei:real;
begin
writeln('->',calcul);
a[1]:=1.0;a[2]:=2.0;
b[1]:=0.0;b[2]:=3.0;
c[1]:=3.0;c[2]:=2.0;
x[1]:=0.2;x[2]:=0.0;
Hypergeometric2F1c(a,b,c,x,rez);
writeln('krit=',rez[1],'+i*',rez[2]);
calcul:=calcul+1;
ef:=send_ef;
ei:=send_ei;
writeln('ei=',ei,'ef=',ef);
writeln('begin fun1');}
{writeln(ki(ei),kf(ef),t);}
{kki:=ki(ei);
kkf:=kf(ef);
writeln(kki,kkf,t);
qq:=2*kki*kkf*(1-t);
qqq:=kki*kki+kkf*kkf-2*t*kki*kkf;
qqqq:=qq/qqq;
writeln(qq,qqq,qqqq);}

{x[1]:=qq/qqq}{2*(1-t)/(ki(ei)/kf(ef)+kf(ef)/ki(ei)-2*t)}{; x[2]:=0.0;}
{writeln('ki=',ki(ei), 'kf=', kf(ef),'   x=',x[1]);}
{writeln('a[2]=',li(ei),'b[2]=',lf(ef));}
{llf:=lf(ef);
lli:=li(ei);}
{a[1]:=1.0;a[2]:=lli;
b[1]:=0.0;b[2]:=-llf;
c[1]:=1.0;c[2]:=0.0;}
{a[1]:=1.0;a[2]:=2.0;
b[1]:=0.0;b[2]:=3.0;
c[1]:=3.0;c[2]:=2.0;
x[1]:=0.2;x[2]:=0.0;
Hypergeometric2F1c(a,b,c,x,rez);
writeln('krit=',rez[1],'+i*',rez[2]);
per:=kki*kki+kkf*kkf;
fun1:=1e-300*(rez[1]*rez[1]+rez[2]*rez[2])/(per-2*kki*kkf*t)/(per-2*kki*kkf*t);}
{writeln('end fun1');}
{end;}

function fun1(t:real):real;
var a,b,c,x,rez:complex;
    lli,llf,kki,kkf,per,ef,ei:real;
begin
writeln('->',calcul);
a[1]:=1.0;a[2]:=2.0;
b[1]:=0.0;b[2]:=3.0;
c[1]:=3.0;c[2]:=2.0;
x[1]:=0.2;x[2]:=0.0;
Hypergeometric2F1c(a,b,c,x,rez);
writeln('krit=',rez[1],'+i*',rez[2]);
calcul:=calcul+1;
ef:=send_ef;
ei:=send_ei;
writeln('ei=',ei,'ef=',ef);
writeln('begin fun1');
{writeln(ki(ei),kf(ef),t);}
kki:=sqrt(2*m*ei);;
kkf:=sqrt(2*m*ef);
writeln(kki,kkf,t);
{qq:=2*kki*kkf*(1-t);
qqq:=kki*kki+kkf*kkf-2*t*kki*kkf;
qqqq:=qq/qqq;
writeln(qq,qqq,qqqq);}

{x[1]:=qq/qqq;}
x[1]:=2*(1-t)/(ki(ei)/kf(ef)+kf(ef)/ki(ei)-2*t);
x[2]:=0.0;
{writeln('ki=',ki(ei), 'kf=', kf(ef),'   x=',x[1]);}
{writeln('a[2]=',li(ei),'b[2]=',lf(ef));}
llf:=(z1+1)*z2*aa*sqrt(m/2/ef);
lli:=z1*z2*aa*sqrt(m/2/ei);
{a[1]:=1.0;a[2]:=lli;
b[1]:=0.0;b[2]:=-llf;
c[1]:=1.0;c[2]:=0.0;}
a[1]:=1.0;a[2]:=2.0;
b[1]:=0.0;b[2]:=3.0;
c[1]:=3.0;c[2]:=2.0;
x[1]:=0.2;x[2]:=0.0;
Hypergeometric2F1c(a,b,c,x,rez);
writeln('krit=',rez[1],'+i*',rez[2]);
per:=kki*kki+kkf*kkf;
fun1:=1e-300*(rez[1]*rez[1]+rez[2]*rez[2])/(per-2*kki*kkf*t)/(per-2*kki*kkf*t);
{writeln('end fun1');}
end;


function fun2(ef:real):real;forward;
function fun3(ei:real):real;forward;

procedure integral(a,b:real;cl:integer;var it:real);
{Процедура расчета интеграла 6-титочечным  методом Гаусса с заданной точностью.
a и b -пределы интегрирования. cl - определяет подинтегральную функцию,int -
значение интеграла}
label 3;
var x,c:array[1..6] of real;
    ing:array[1..2] of real;
    i,k,n,j:integer;
    pa,pb,pf,h:real;

begin
 c[1]:=0.085662246;
 c[2]:=0.180380786;
 c[3]:=0.233956967;
 c[6]:=c[1];
 c[5]:=c[2];
 c[4]:=c[3];
 x[1]:=-0.932469514;
 x[2]:=-0.661209386;
 x[3]:=-0.238619186;
 x[6]:=-x[1];
 x[5]:=-x[2];
 x[4]:=-x[3];

 n:=1;
 3: ing[1]:=0;
    ing[2]:=0;
    for k:=1 to 1 do
 begin
     h:=(b-a)/n;
     for i:=1 to n do
    begin
      pa:=a+h*(i-1);
      pb:=a+h*i;
      for j:=1 to 6 do
      begin
      if cl=1 then
      pf:=(pb-pa)*fun1((pb+pa)/2+(pb-pa)*x[j]/2)/2;
      if cl=2 then
      pf:=(pb-pa)*fun2((pb+pa)/2+(pb-pa)*x[j]/2)/2;
      if cl=3 then
      pf:=(pb-pa)*fun3((pb+pa)/2+(pb-pa)*x[j]/2)/2;

      ing[k]:=ing[k]+2*c[j]*pf;
       end;
     end;
    n:=n*2;
 end;
{if (abs(ing[2]-ing[1]))>toch then goto 3;}
it:=ing[1];
end;

function sig1(ef:real):real;
var igf:real;
begin
send_ef:=ef;
{writeln('begin sig1');}
integral(-1.0,1.0,1,igf);
sig1:=igf;
{writeln('end sig1');}
end;

function fun2(ef:real):real;
var e,c:extended;
      ei:real;
begin
   ei:=send_ei;{writeln('begin fun2');}
   e:=ei-d-df-1;
   c:=1e300*sqrt(e-ef)*(e-ef)*(e-ef)*(e-ef)/(ki(ei)*ki(ei)-kf(ef)*kf(ef))/(ki(ei)*ki(ei)-kf(ef)*kf(ef));
   fun2:=c/(exp(2*pi*lf(ef))-1)*sig1(ef);
   {writeln('end fun2');}
end;

function sigb(ei:real):real;
var igr:real;
begin
{writeln('begin sigb');}
   integral(0.0,ei-d-df-1,2,igr);
   sigb:=z2*z2*256*sqrt(2)*aa*aa*aa*aa*gv*gv*m*m*m*m*m*z1*(z1+1)*z2*z2/(105*pi*ei)/(1-exp(-2*pi*li(ei)))*igr;
end;

function fun3(ei:real):real;
begin
   writeln('begin fun3');
   send_ei:=ei;
   fun3:=exp(-ei/tt)*ei*sigb(ei);
end;

function sigv(d:real):real;
var igt:real;
begin
{writeln('begin sigb');}
   integral(d+df+1,100,3,igt);
   sigv:=sqrt(8/pi/m/tt/tt/tt)*igt*0.09747;
end;

begin

calcul:=1;

clrscr;
77: writeln;
writeln('Введите число точек, для которых будет производиться расчет:');
write('n=');
readln(n);
assign(f,'d:\u\pascal\t\sigv\nsv1.dat');
append(f);
close(f);
append(f);
for i:=1 to n do begin
writeln('Введите значения:');
write('d=');
readln(d);
write(f,'d=',d:10);
d:=d/mev ;
write('df=');
readln(df);
write(f,' df=',df:10);
df:=df/mev;
write('z1=');
readln(z1);
writeln(f,' z1=',z1:10);
write('z2=');
readln(z2);
writeln(f,' z2=',z2:10);
write('a1=');
readln(a1);
writeln(f,' a1=',a1:10);
write('a2=');
readln(a2);
writeln(f,' a2=',a2:10);
write('tt=');
readln(tt);
writeln(f,' tt=',tt:10);
tt:=tt*1e-10/1.16/mev;

mn:=1.835E3;
gv:=2.9899e-12;
aa:=1/137.03;
toch:=1e3;
m:=mn*a1*a2/(a1+a2);
writeln(a1*a2/(a1+a2));


{for id:=0 to 0 do
 begin
  ei:=100.0+0.5*id;
  ei:=ei/mev;
  ki:=sqrt(2*m*ei);
  li:=z1*z2*aa*sqrt(m/2/ei);
 { writeln('ei=',ei*mev,'     ','sigb=',sigb(ei)*1.4915e-21);}
  writeln(f,'tem=',tt*1.16*mev,'     ','sigv=',sigv(d)*44.722E-12:14);
  writeln('OK');

writeln;
end;
{end;}
close(f);
writeln('Необходимо ли Вам продолжить расчеты?');
an:=readkey;
if an='y' then goto 77;
end.




