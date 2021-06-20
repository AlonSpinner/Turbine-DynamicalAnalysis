%% initalize
mili=1e-3; micro=1e-6;
E=210*1e9; %Pa
L=65; %m
m=90*1e3; %kg
r_out=3; %m
r_in=2.5; %m
c1=22*1e4; %Ns/m, Fd=-c1*du/dt

%% A1
I=pi*(r_out^4-r_in^4)/(4); %m^4
k=(3*E*I)/(L^3); %N/m

%% A2 - no code needed.
%motion eq: mu''+c1u'+ku=0
%' represents time derivitive

%% A3
wn=sqrt(k/m); %rad/s
zeta=c1/(2*sqrt(k*m)); %unitless
wd=wn*sqrt(1-zeta^2); %rad/s. free osciliation radial frequency
wr=wn*sqrt(1-2*zeta^2); %rad/s. maximal output/input 4 steady state radial frequency

%% A4
fn=wn/(2*pi);
t=linspace(0,3/fn,100); %s, we want to show 3 machzorim
beta=linspace(0,2*wn,100);
[T,Beta]=meshgrid(t,beta); %mesh matrix form coordinates
% u_InShock=(1040*Beta.*sqrt(Beta)./(wn*pi*(Beta.^2-wn.^2))).*(wn*sin(Beta.*T)-Beta.*sin(wn.*T));
u=(1040*Beta.*sqrt(Beta)./(wn*pi*(Beta.^2-wn.^2))).*cos((pi*wn)./(2*Beta)).*sin(wn.*T-(pi*wn)./(2.*Beta));

%Plot u(beta,t)
srffig=figure;
srfAx=axes(srffig);
surf(srfAx,T,Beta,u)
xlabel(srfAx,'t [s]'); 
ylabel(srfAx,'\beta [rad/s]'); 
zlabel(srfAx,'u [m]');

%plot max(u(beta,t),beta) and wn
Maxu=max(u,[],2);
Srspfig=figure;
SrspAx=axes(Srspfig); grid(SrspAx,'on'); hold(SrspAx,'on');
plot(SrspAx,beta,Maxu,'linewidth',2)
plot(SrspAx,[wn,wn],SrspAx.YLim,'--g','linewidth',2)
xlabel(SrspAx, '\beta [rad/s]'); 
ylabel(SrspAx, 'Displacement [m]');
legend('max(u(t,\beta))','w_n','location','best');

%% PA5
%Eq: mu''+c1u'+ku=F;   F=F0cos(wt)
H=tf(1,[m,c1,k]);
BodeFig=figure;
BodeAx=axes(BodeFig);
bode(BodeAx,H);
grid(BodeAx,'on');

%% A6
A=3.5*1e6; %N
B=2.5*1e6; %N;
alpha=0.1; %rad/m
phi=pi/3; %rad

aUinf=linspace(0.1*wn,2*wn,100); %rad/s
Uinf=aUinf/alpha; %m/s
t=linspace(0,2*pi/(0.1*wn),100); %build time vector according to slowlest frequency
[T,AUINF]=meshgrid(t,aUinf);
r=aUinf/wn; %Normalize
H=@(r) 1./(k*(1-r.^2+2*1i*zeta*r)); %normalized transfer function form
M1=(abs(H(r)))'; Psy1=(angle(H(r)))'; %hegber and phasa of alpha*Uinf
M3=(abs(H(3*r)))'; Psy3=(angle(H(3*r)))'; %hegber and phasa of 3*alpha*Uinf
U1=A*M1.*cos(AUINF.*T+Psy1); %m %AUINF changed across rows, and so Psy1 is column vector
U2=B*M3.*cos(AUINF.*T+Psy3+phi); %m
U=U1+U2;
MaxU=max(U,[],2);
MaxSigma=k*L*r_out*MaxU/I; %Pa
Sigma_allowed=108; %Mpa

%Plot Stress vs windspeed
sigmafig=figure;
sigmaAx=axes(sigmafig); hold(sigmaAx,'on'); grid(sigmaAx,'on');
plot(sigmaAx,Uinf,MaxSigma*micro,'linewidth',2);
plot(sigmaAx,sigmaAx.XLim,[Sigma_allowed,Sigma_allowed],'r','linewidth',2);
plot(sigmaAx,[wn/alpha,wn/alpha],sigmaAx.YLim,'--g','linewidth',1)
plot(sigmaAx,[wn/(3*alpha),wn/(3*alpha)],sigmaAx.YLim,'--m','linewidth',1)
xlabel(sigmaAx,'U_\infty [m/s]'); 
ylabel(sigmaAx,'\sigma [Mpa]'); 
legend(sigmaAx,'Maximal stress in beam','Allowed stress','w_n/\alpha','(1/3)w_n/\alpha');

%% B1
%nothing to compute. motion equation 

%% B2
%find potenial(U) and its derivative with respect to time (dU)
%find equilibrium points from dU=0, classify them from U

%initalize new parameters
ks=4*1e13; %N/m
c5=1e3; %N*s^5/m^5
g=@(y) (1/25)*(0.5*y.^4-0.25*y.^2);
dg=@(y) (1/25)*(2*y.^3-0.5*y);

%build u,U,dU
u=linspace(-0.8,0.8,1000); %m
U=@(u) 0.5*k*u.^2+ks*(g(u)).^2; %J
dU=@(u) k*u+2*ks.*g(u).*dg(u);

%find equilibriums
ForceSpline=csapi(u,dU(u)); %Cubic Spline interpolation - builds piece wise 3d degree polys. 
%can be plotted with fnplt(spline)
ZeroInt=fnzeros(ForceSpline,[min(u),max(u)]); %find zeros of the spline in a given interval
%Z is a 2xm numeric vector array: each column 1:m describes the interval over which the spline is zero.
%if the spline reaches zero, but does not cross the zero line the root may
%not reveal itself in ZeroInt.
ZeroPoint=abs(ZeroInt(2,:)-ZeroInt(1,:))<eps; %bool array. 1: only the zero intervals which are of 0 length == a root.
Roots=ZeroInt(1,ZeroPoint); %take the first row (arbitrary) from ZeroInt which complies with ZeroPoint

%do some plotting
potentfig=figure;
potentAx=axes(potentfig); hold(potentAx,'on'); grid(potentAx,'on');
plot(potentAx,u,U(u)*micro,'linewidth',2);
scatter(potentAx,Roots,U(Roots)*micro,50,'filled');
plot([0.5,0.5],potentAx.YLim,'--m')
plot([-0.5,-0.5],potentAx.YLim,'--m')
xlabel(potentAx,'u [m]'); ylabel(potentAx,'Potential Energy [MJ]');
legend(potentAx,'Potential Energy','Equilibrium points','peak of g(u)',...
    'location','north');

%% B3
%equation of motion around equlibrium:
%m*ddu+Ceq*du+Keq*u=F

check=@(u_eq) (56*u_eq^6-30*u_eq^4+3*u_eq^2)/8;
Keq=@(u_eq) k+2/(25^2)*ks*(56*u_eq^6-30*u_eq^4+3*u_eq^2)/8;
Ceq=@(u_eq) c1; %du in equilbrium is 0

%for u_eq1=0 the results are the same as in A2

%for u_eq2=0.7037
wn_eq2=sqrt(Keq(0.7037)/m); %rad/s
zeta_eq2=Ceq(0.7037)/(2*sqrt(Keq(0.7037)*m)); %unitless
wd_eq2=wn_eq2*sqrt(1-zeta_eq2^2); %rad/s. free osciliation radial frequency
wr_eq2=wn_eq2*sqrt(1-2*zeta_eq2^2); %rad/s. maximal output/input 4 steady state radial frequency
%% B4
%should I use ODE45 to solve linear equation, or non linear equation?
A=3.5*1e6; %N
B=2.5*1e6; %N;
alpha=0.1; %rad/m
phi=pi/3; %rad

%first stable equalbrium u0=0. second stable equilibrium is it self in
%"danger zone" stress wise. no point in checking
u0=[0,0];
aUinf=linspace(0.1*wn,2*wn,100); %rad/s
T=2*pi/(0.1*wn); %slowest frquency
tspan=linspace(0,20*T,10000); %tons of cycles
[~,du]=arrayfun(@(aUinf) ode45(@(t,u) odefunB6(t,u,m,c1,c5,k,ks,g,dg,A,B,phi,aUinf),tspan,u0),...
    aUinf,'un',0); %du = [u,dudt] output
du=cellfun(@(x) x(:,1),du,'un',0); %obtain all odd columns
u=cell2mat(du)-u0(1);
u_ss=u(9000:10000,:);
MaxU=(max(abs(u_ss),[],1))';
Uinf=aUinf/alpha; %m/s
MaxSigma=k*L*r_out*MaxU/I; %Pa
Sigma_allowed=108; %Mpa

%Plot Stress vs  windspeed
sigmafig=figure;
sigmaAx=axes(sigmafig); hold(sigmaAx,'on'); grid(sigmaAx,'on');
plot(sigmaAx,Uinf,MaxSigma*micro,'linewidth',2);
plot(sigmaAx,sigmaAx.XLim,[Sigma_allowed,Sigma_allowed],'r','linewidth',2);
xlabel(sigmaAx,'U_\infty [m/s]'); 
ylabel(sigmaAx,'\sigma [Mpa]'); 
legend(sigmaAx,'Maximal stress in beam','Allowed stress');

function dudt=odefunB6(t,u,m,c1,c5,k,ks,g,dg,A,B,phi,aUinf)
%equation from linearization around stable equilibrium
%syntax:
%u(1)=u; u(2)=dudt
%dudt(1)=du/dt; dudt(2)=d2u/dt2

dudt=zeros(2,1);
dudt(1)=u(2);
dudt(2)=(1/m)*(-c1*u(2)-c5*(u(2))^5-k*u(1)-2*ks*g(u(1))*dg(u(1))+A*cos(aUinf*t)+B*cos(3*aUinf*t+phi));
end