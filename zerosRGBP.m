function [zzero,ier]=zerosRGBP(n,a)
%---------------------------------------------
% Computation of the complex zeros of the
% reverse generalized Bessel polynomials
%------------------------------------------------------
% Inputs:
%     n,a: parameters of the function. 
%          n values are expected to be moderate/large. 
% Outputs:
%     zzero: array with the computed zeros.
%     ier:   error flag indicating the success or issues 
%            in computation:
%               0: Computation successful
%               1: Expected loss of accuracy for one or 
%                  more of the computed zeros. 
%----------------------------------------------------------------
%  Matlab implementation of the algorithm described 
%  in the paper
%  "Uniform asymptotic approximation and numerical
%   evaluation of the reverse generalized Bessel
%   polynomial zeros", 
%  by T.M. Dunster, A. Gil, D. Ruiz-Antolin, J. Segura
%---------------------------------------------------------------- 
epsilon=1e-12;
ier=0;
niter=100;
cc=2*(-floor(n/2)+floor(n)/2);
mm=floor((n-cc)/2)-(1-cc);
m=mm+1;
zinit=zasympRGBP(a,n,m);
zc(1)=zinit;
j=1;
istop=0;
while (istop==0)  
  h=Hp(a,n,zc(j));
  z0=zc(j);
  z=z0+h;
  y0=0;
  y1=1;
  delta=1+epsilon;
  it=0;
  while (delta>epsilon && it<niter)
    y=z; 
    [z,hn,s,sd]=Tptayl(a,n,z0,h,y0,y1);
    y0=s;
    y1=sd;
    h=hn;
    z0=z;
    zz=z0+h;
    delta=abs(zz-y)/abs(y);
    it=it+1;
  end
  if it==niter
    ier=1;
  end
  if isnan(zz)==0
    j=j+1;
    zc(j)=zz;
  else  
    istop=1;  
  end  
  if j==m
    istop=1;
  end
end
kl=0;
for i=1:j
  if real(zc(i))<0
    kl=kl+1;
    zzero(kl)=zc(i);
  end  
end
zzero=zc;
end

function h=Hp(a,n,z)   
Omega=-1+(2-a)/z-(n+a/2)*(n+a/2-1)/(z*z);
coef=sqrt(Omega);
h=pi/coef;
end

function [z,hn,s,sd]=Tptayl(a,n,zz,h,y0,y1)
[s,sd]=step(a,n,zz,h,y0,y1);
H=s/sd;
z=zz+h;
Omega=-1+(2-a)/z-(n+a/2)*(n+a/2-1)/(z*z);
coef=sqrt(Omega);
hn=-1/coef*atan(coef*H);
end

function [s,sd]=step(a,n,z,h,y0,y1)
% Auxiliary routine for taylor
imax=30;       
onesix=1/6; 
xin=z;
xin2=xin*xin;
xin3=xin2*xin;
h2=h*h;
h2d2=h2*0.5;
h3=h2*h;
xifac=1-(2-a)/xin+(2*n+a)*(2*n+a-2)/(4*xin2);
xifacp=(2-a)/xin2-2*(n+1/2*a)*(n+1/2*a-1)/xin3;
y2=xifac*y0;
y3=xifac*y1+xifacp*y0;
s=y0+y1*h+y2*h2d2+y3*h3*onesix;
sd=y1+y2*h+y3*h2d2;
lt=h3*onesix;
ltd=h2*0.5;
for i=0:imax
  di=i;
  k=i+2;
  coun=di+4.0;
  cound=di+3.0;
  y4=-(2*k*xin*y3+(-xin2+xin*(2-a)-(2*n+a)*(2*n+a-2)/4+k*(k-1))*y2+...
       k*(-2*xin+2-a)*y1-k*(k-1)*y0)/xin2; 
  lt=lt*h/coun;
  ltd=ltd*h/cound;
  s=s+lt*y4;
  sd=sd+ltd*y4;
  y0=y1;
  y1=y2;
  y2=y3;
  y3=y4;
end
end

function zinit=zasympRGBP(a,n,m)
u=n+0.5;
akz=[-2.33810741045976703848919725245;
    -4.08794944413097061663698870146;
    -5.52055982809555105912985551293;
    -6.78670809007175899878024638450;
    -7.94413358712085312313828055580;
    -9.02265085334098038015819083988;
    -10.0401743415580859305945567374;
    -11.0085243037332628932354396496;
    -11.9360155632362625170063649029;
    -12.8287767528657572004067294072;
    -13.6914890352107179282956967795;
    -14.5278299517753349820739814430;
    -15.3407551359779968571462085135;
    -16.1326851569457714393459804472;
    -16.9056339974299426270352387706;
    -17.6613001056970575092536503040;
    -18.4011325992071154158613979295;
    -19.1263804742469521441241486897;
    -19.8381298917214997009475636160;
    -20.5373329076775663599826814113];
if m>20
  t=3.0*pi/8*(4.0*m-1.0);
  t2=t*t;
  t4=t2*t2;
  t6=t4*t2;
  t8=t6*t2;
  t10=t8*t2;
  akm=-t^(2.0/3.0)*(1.0+5/(48*t2)...
       -5.0/(36*t4)+77125.0/(82944*t6)...
       -108056875.0/(6967296*t8)+162375596875.0/(334430208.0*t10));
else
  akm=akz(m);  
end
zeta=akm/u^(2/3);
xi=-2/(3*u)*1i*abs(akm)^(3/2);
rho=1/zeta;
alpha=(a-2)/u;
sigma=sqrt(1+alpha);
sigma2=sigma*sigma;
sigma3=sigma2*sigma;
sigma4=sigma2*sigma2;
sigma5=sigma3*sigma2;
sigma6=sigma4*sigma2;
sigma7=sigma5*sigma2;
sigma8=sigma6*sigma2;
sigma9=sigma8*sigma;
sigma10=sigma8*sigma2;
z=solveq316(xi,alpha);
Z=-1/2*sqrt(alpha^2 + 4*alpha*z + 4*z^2 + 4*alpha + 4);
C= (z+alpha/2)/Z;
C2=C*C;
C3=C2*C;
C4=C3*C;
C5=C4*C;
C6=C5*C;
C7=C6*C;
C8=C7*C;
C9=C8*C;
C10=C9*C;
C11=C10*C;
C12=C11*C;
C13=C12*C;
C14=C13*C;
C15=C14*C;
S= sigma/Z;
xi1= sigma/(z*S);
xi2= -sigma/(z^2*S)+C/z;
rho1=-3/2*rho^4*xi*sigma/(z*S);
rho2=9*rho^7*xi^2*sigma2/(z^2*S^2)-(3/2)*rho^4*sigma2/(z^2*S^2)+...
     (3/2)*rho^4*xi*sigma/(z^2*S)-(3/2)*rho^4*xi*C/z;
zeta1= (2/3)*zeta*sigma/(z*S*xi);
zeta2= -(2/9)*zeta*sigma2/(z^2*S^2*xi^2)-(2/3)*zeta*sigma/(z^2*S*xi)+...
         (2/3)*zeta*C/(z*xi);
zeta3= (2/3)*zeta*(C6*xi^2*z^2-3*C4*xi^2*z^2-...
       2*sigma*xi*z*(xi*S+1/2*sigma)*C3+xi^2*(2*sigma2+...
       3*z^2)*C2+2*sigma*xi*z*(xi*S+1/2*sigma)*C-...
       sigma3*S*xi+(-2*sigma2-z^2)*xi^2-...
       4/9*sigma4)/S/z^3/xi^3/sigma/(C2-1);
alpha2=alpha*alpha;
alpha3=alpha2*alpha;
alpha4=alpha3*alpha;
alpha5=alpha4*alpha;
alpha6=alpha5*alpha;
alphap1=1+alpha;
alphap12=alphap1*alphap1;
alphap13=alphap12*alphap1;
alphap15=alphap13*alphap1;
alphap17=alphap15*alphap1;
d(1)=-1/48*alpha/alphap1;
d(3)=(7/5760)*alpha*(alpha2+3*alpha+3)/alphap13;
d(5)=-31/80640*alpha*(alpha4+5*alpha3+10*alpha2+10*alpha+5)/alphap15;
d(7)=(127/430080)*alpha*(alpha6+7*alpha5+21*alpha4+35*alpha3+...
     35*alpha2+21*alpha+7)/alphap17;
E(1)=(1/48)*((5*sigma2-5)*C3+10*S*C2*sigma+...
    (-6*sigma2+6)*C+sigma2-4*S*sigma-1)/sigma2;
E(3)=(1/46080)*((5525*sigma6-82875*sigma4+82875*sigma2-...
      5525)*C^9+33150*S*(sigma2-3)*sigma*(sigma2-1/3)*C8+...
      (-19890*sigma6+258570*sigma4-258570*sigma2+19890)*C7-...
      92820*S*(sigma4-64/21*sigma2+1)*sigma*C6+...
      (26541*sigma6-288855*sigma4+288855*sigma2-26541)*C5+...
      88470*S*sigma*(sigma4-40282/14745*sigma2+1)*C4+...
      (-15540*sigma6+133500*sigma4-133500*sigma2+15540)*C3-...
      31080*S*sigma*(sigma4-9218/3885*sigma2+1)*C2+(3420*sigma6-...
      20340*sigma4+20340*sigma2-3420)*C+(2280*sigma5-4432*sigma3+...
      2280*sigma)*S-56*sigma6+56)/sigma6;
E(5)=1/10321920*((8696625*sigma10-391348125*sigma8+...
     1826291250*sigma6-1826291250*sigma4+391348125*sigma2-...
     8696625)*C15+86966250*S*sigma*(sigma4-2*sigma2+...
     1/5)*(sigma4-10*sigma2+5)*C14-...
     52179750*(sigma-1)*(sigma+1)*(sigma8-40*sigma6+...
     142*sigma4-40*sigma2+1)*C13-452224500*(sigma8-...
     144/13*sigma6+294/13*sigma4-144/13*sigma2+1)*S*sigma*C12+...
     (132029730*sigma10-4874921730*sigma8+20487621420*sigma6-...
     20487621420*sigma4+4874921730*sigma2-132029730)*C11+...
     968218020*(sigma8-7783848/768427*sigma6+...
     15411826/768427*sigma4-7783848/768427*sigma2+...
     1)*S*sigma*C10+(-181835360*sigma10+5954923800*sigma8-...
     23554911240*sigma6+23554911240*sigma4-5954923800*sigma2+...
     181835360)*C^9-1091012160*(sigma8-17829949/1948236*sigma6+...
     5688077/324706*sigma4-17829949/1948236*sigma2+1)*S*sigma*C8+...
     (146258505*sigma10-4160688525*sigma8+15372246330*sigma6-...
     15372246330*sigma4+4160688525*sigma2-146258505)*C7+...
     682539690*S*sigma*(sigma8-79312556/9750567*sigma6+...
     341039042/22751323*sigma4-79312556/9750567*sigma2+...
     1)*C6+(-67995018*sigma10+1629504450*sigma8-...
     5564884500*sigma6+5564884500*sigma4-1629504450*sigma2+...
     67995018)*C5-226650060*(sigma8-3809876/539643*sigma6+...
     78411122/6295835*sigma4-3809876/539643*sigma2+1)*S*sigma*C4+...
     (16611420*sigma10-318811500*sigma8+990816120*sigma6-...
     990816120*sigma4+318811500*sigma2-16611420)*C3+...
     33222840*S*sigma*(sigma8-699548/118653*sigma6+...
     13661566/1384285*sigma4-699548/118653*sigma2+1)*C2+...
     (-1590120*sigma10+21971880*sigma8-60464880*sigma6+...
     60464880*sigma4-21971880*sigma2+1590120)*C+(-1060080*sigma9+...
     4820480*sigma7-7528992*sigma5+4820480*sigma3-...
     1060080*sigma)*S+3968*sigma10-3968)/sigma10;
Ep(1)=5/16*S^2*((C2-2/5)*(sigma-1)*(sigma+1)*S-2*(C2-...
      4/5)*C*sigma)/sigma3;
Epp(1)=-25/16*S^3*(-2*C4*sigma+S*(sigma2-1)*C3+...
       54/25*C2*sigma+S*(-16/25*sigma2+16/25)*C-8/25*sigma)/sigma4;
Ep(3)=-1105/1024*S^4*((-6*sigma5+20*sigma3-6*sigma)*C7+...
       S*(sigma6-15*sigma4+15*sigma2-1)*C6+(62/5*sigma+...
       62/5*sigma5-188/5*sigma3)*C5-9/5*S*(sigma4-...
       98/9*sigma2+1)*(sigma+1)*(sigma-1)*C4-8504/1105*(sigma4-...
       2878/1063*sigma2+1)*sigma*C3+192/221*S*(sigma4-...
       39/5*sigma2+1)*(sigma+1)*(sigma-1)*C2+1432/1105*sigma*...
       (sigma4-422/179*sigma2+1)*C-76/1105*S*(sigma4-94/19*sigma2+...
       1)*(sigma+1)*(sigma-1))/sigma7;
Epp(3)=12155/1024*S^5*((-6*sigma5+20*sigma3-6*sigma)*C8+...
       S*(sigma6-15*sigma4+15*sigma2-1)*C7+768/55*(sigma4-...
       299/96*sigma2+1)*sigma*C6-111/55*S*(sigma4-434/37*sigma2+...
       1)*(sigma+1)*(sigma-1)*C5-128038/12155*(sigma4-...
       184454/64019*sigma2+1)*sigma*C4+14676/12155*(sigma4-...
       34762/3669*sigma2+1)*S*(sigma+1)*(sigma-1)*C3+...
       32672/12155*(sigma4-2686/1021*sigma2+1)*sigma*C2-...
       460/2431*(sigma4-4214/575*sigma2+1)*S*(sigma+1)*(sigma-1)*C-...
       1432/12155*sigma*(sigma4-422/179*sigma2+1))/sigma8;
Y(1)= 3/2*xi*(E(1)+d(1))*rho^2-5/48*rho^2;
Y(2)= -1/4*Y(1)^2*rho+3/2*xi*(E(3)+d(3))*rho^2 +...
      5/32*Y(1)*rho^3-1105/9216*rho^5;
Y(3)= 1/24*Y(1)^3*rho^2-25/128*Y(1)^2*rho^4-...
     1/2*Y(1)*Y(2)*rho+1105/2048*Y(1)*rho^6+...
     3/2*xi*(E(5)+d(5))*rho^2+5/32*Y(2)*rho^3-82825/98304*rho^8;
Yp(1)= 1/24*rho*(36*rho*xi*Ep(1) + 36*rho*xi1*(E(1)+d(1))+...
       72*rho1*xi*(E(1)+d(1))-5*rho1);
Ypp(1)= 1/24*(36*xi*Epp(1)+72*xi1*Ep(1)+36*xi2*(E(1)+d(1)))*rho^2+...
       1/24*((144*rho1*xi1+72*rho2*xi)*(E(1)+d(1))+...
       144*xi*rho1*Ep(1)-5*rho2)*rho + 3*xi*(E(1)+d(1))*rho1^2-...
       5/24*rho1^2;
Yp(2)= -1/2*Y(1)*rho*Yp(1)-1/4*Y(1)^2*rho1+3/2*xi1*(E(3)+d(3))*rho^2+...
        3/2*xi*Ep(3)*rho^2+3*xi*(E(3)+d(3))*rho*rho1+...
        5/32*Yp(1)*rho^3+15/32*Y(1)*rho^2*rho1-5525/9216*rho^4*rho1;
tau(1)=z;
tau(2)= -Y(1)/zeta1;
tau(3)= -1/2*(zeta2*tau(2)^2 + 2*Yp(1)*tau(2) + 2*Y(2))/zeta1;
tau(4)= -1/6*(tau(2)^3*zeta3+3*tau(2)^2*Ypp(1)+...
       6*zeta2*tau(2)*tau(3)+6*Yp(1)*tau(3)+6*Yp(2)*tau(2)+6*Y(3))/zeta1;
zapp=0;
for k=0:3
  zapp=zapp+tau(k+1)/u^(2*k);
end
zinit=u*zapp;
end

function zsol=solveq316(xi,alpha)
epsi=1e-14;
niter=100;
w0=-0.1-0.1*1i;
er=1;
iter=0;
while er>epsi && iter<niter
  w=w0-feq316(xi,alpha,w0)/feq316p(alpha,w0);
  er=abs(1-w/w0);
  w0=w;
  iter=iter+1;
end
zsol=-0.5+w0;
end  

function f=feq316(xi,alpha,w)
beta=(1/2)*log(1+alpha)+(2+(1/2)*alpha)*log(2)-...
     (1/2)*1i*pi;     
f=-(1/2)*(4*(-.5+w)^2+4*alpha*(-.5+w)+...
 alpha^2+4*alpha+4)^(1/2)+(1+(1/2)*alpha)*...
 log((-.5+w)/(alpha^2+2*(-(1/2)*(4*(-.5+w)^2+...
  4*alpha*(-.5+w)+alpha^2+4*alpha+4)^(1/2)+1.5+w)*alpha-...
  2*(4*(-.5+w)^2+4*alpha*(-.5+w)+alpha^2+4*alpha+4)^(1/2)+4))+...
  (1/2)*alpha*log(1.0-2*w-alpha+(4*(-.5+w)^2+...
  4*alpha*(-.5+w)+alpha^2+4*alpha+4)^(1/2))+beta-xi;
end

function fp=feq316p(alpha,w)
fp=-(5.00-4.0*w+4.*w^2+2.0*alpha+...
  4.*alpha*w+alpha^2)^(1/2)/(-1.+2.*w);
end
