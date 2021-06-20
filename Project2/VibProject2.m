%% A-initalize
mili=1e-3; micro=1e-6;
E=210e9; %Pa
L=65; %m
l=2; %m
m=90e3; %kg
r_out=3; %m
r_in=2.5; %m
c=22e4; %Ns/m 
I=pi*(r_out^4-r_in^4)/(4); %m^4
%% A-Q1
%Mass and damping matrices
Mmat=m*eye(2);
Cmat=c*eye(2);

%Stiffness Matrix
syms x sL sl u v sEI real
p=x.^(0:3); %build [1,4] vector polynomial [1,x,x^2,x^3]
%Solve deflection function w
A=[subs(p,x,0); subs(p,x,sL);
    subs(diff(p,x),x,0); subs(diff(p,x),x,sL)];
b=[0;u;0;-atan2(v,sl)];
a=A\b;
w=p*a;

%Calculate stiffness matrix K
ddw=diff(w,x,2);
V=0.5*sEI*int(ddw^2,x,0,sL);
K = simplify([  diff(V,v,2)  diff(diff(V,u),v) 
        diff(diff(V,v),u)   diff(V,u,2)]);
%% A-Q2
u0=0; v0=0;
K0 = double(subs(K,[sL,sl,v,u,sEI],[L,l,v0,u0,E*I]));
[PHI,wr_2]=eig(K0,Mmat); %solves Av=lambda*Bv ~ K-(wn^2)*M
wr=sqrt(diag(wr_2));
%% A-Q3
%noramlized mass
Gama=PHI'*Cmat*PHI;
Gama=0.5*(Gama+Gama'); %enfroce sym
Lambda=PHI'*K0*PHI;
Lambda=0.5*(Lambda+Lambda'); %enforce sym
%% A-Q5
%Initial conditions 
q0=[3e-4,-66.6e-4,0,0;
    0,0,-33.3,-1.5;
    -15.9e-4,-17.4e-4,0,0]';
InvPHI=PHI'*Mmat; %faster than inv(phi)
eta0=[InvPHI*q0(1:2,:);InvPHI*q0(3:4,:)];
%define state vector z as [eta1,eta2,deta1,deta2]
%so dz=[deta1,ddeta1,deta2,ddeta2]
Ad=[zeros(2),eye(2);-Lambda,-Gama];
Cd=[eye(2),zeros(2)];
sys=ss(Ad,[],Cd,[]);

%Build time vector to sum upon
Cycles=20;
T1=2*pi/wr(1);
T2=2*pi/wr(2);
Fs=10*(1/T2);
t=0:1/Fs:(Cycles*T1);
[etaInt1,time1]=initial(sys,eta0(:,1),t);
[etaInt2,time2]=initial(sys,eta0(:,2),t);
[etaInt3,time3]=initial(sys,eta0(:,3),t);

%plots
eta0fig=figure;
eta11Ax=subplot(3,2,1,'parent',eta0fig);
eta12Ax=subplot(3,2,3,'parent',eta0fig);
eta13Ax=subplot(3,2,5,'parent',eta0fig);
eta21Ax=subplot(3,2,2,'parent',eta0fig);
eta22Ax=subplot(3,2,4,'parent',eta0fig);
eta23Ax=subplot(3,2,6,'parent',eta0fig);
hold(eta11Ax,'on'); hold(eta21Ax,'on'); 
hold(eta12Ax,'on'); hold(eta22Ax,'on');
hold(eta13Ax,'on'); hold(eta23Ax,'on');
grid(eta11Ax,'on'); grid(eta21Ax,'on'); 
grid(eta12Ax,'on'); grid(eta22Ax,'on');
grid(eta13Ax,'on'); grid(eta23Ax,'on');
xlabel(eta11Ax,'Time [s]'); ylabel(eta11Ax,'\eta_1 [m]');
xlabel(eta12Ax,'Time [s]'); ylabel(eta12Ax,'\eta_1 [m]');
xlabel(eta13Ax,'Time [s]'); ylabel(eta13Ax,'\eta_1 [m]');
xlabel(eta21Ax,'Time [s]'); ylabel(eta21Ax,'\eta_2 [m]');
xlabel(eta22Ax,'Time [s]'); ylabel(eta22Ax,'\eta_2 [m]');
xlabel(eta23Ax,'Time [s]'); ylabel(eta23Ax,'\eta_2 [m]');
lineColor=lines(3);
plot(eta11Ax,time1,etaInt1(:,1),'color',lineColor(1,:)); plot(eta21Ax,time1,etaInt1(:,2),'color',lineColor(1,:)); 
plot(eta12Ax,time2,etaInt2(:,1),'color',lineColor(2,:)); plot(eta22Ax,time2,etaInt2(:,2),'color',lineColor(2,:)); 
plot(eta13Ax,time3,etaInt3(:,1),'color',lineColor(3,:)); plot(eta23Ax,time3,etaInt3(:,2),'color',lineColor(3,:)); 
legend(eta11Ax,'First Intial condition'); legend(eta21Ax,'First Intial condition');
legend(eta12Ax,'Second Intial condition'); legend(eta22Ax,'Second Intial condition'); 
legend(eta13Ax,'Third Intial condition'); legend(eta23Ax,'Third Intial condition');
%% A-Q6
%P=B*cos(3*aUinf*t+Psy)
%F=Acos(aUinf*t)
%Q=[P,F]
%q=sum_r(phi*H_eta*phi'*Q)

A=1.9e6; %N
B=0.8e6; %N;
alpha=10; %rad/m
psi=pi/3; %rad

aUinf=linspace(0.1*min(wr),2*max(wr),10000); %rad/s
Uinf=aUinf/alpha; %m/s
%Some added functions for simulation
H=@(w) diag(1./diag((-w^2)*eye(2)+(1i*w)*Gama+Lambda)); %2x2 matrix, %normalized mass
% from subs(ddw,[x,sl,sL],[0,l,L]); we obtained bending moment as a
% function of v,u
Momentx0=@(v,u) -E*I*((6*u)/4225+(2*atan2(v,2))/65); %Nm
sigma=@(Mmax) abs(Mmax)*r_out/I; %Pa
MaxSigma=zeros(size(aUinf));

for kk=1:length(aUinf) %loop around and build MaxSigma array
    %Generelized forces Q=[P;F] P=@(w,t) B*exp(3i*w*t+1i*psi); F=@(w,t) A*exp(1i*w*t);
    %super position of inputs for P and F
    %eta1,2 are complex vectors of Amplitude*exp(i*phase)
    eta1=(H(3*aUinf(kk))*PHI'*[B*exp(psi*1i);0]); %for P input
    eta2=(H(aUinf(kk))*PHI'*[0;A]); %for F input
    %Build time vector to sum upon
    Cycles=5;
    T=2*pi/(aUinf(kk));
    Fs=10*(1/T);
    t=0:1/Fs:(Cycles*T);
    %Calculate generalized co-ordinates q 
    q=PHI*real((eta1*exp(3i*aUinf(kk)*t)+eta2*exp(1i*aUinf(kk)*t)));
    MomentMax=arrayfun(@(v,u) Momentx0(v,u),q(1,:),q(2,:));
    MaxSigma(kk)=max(sigma(MomentMax));
    disp([num2str(kk),'/',num2str(length(aUinf))]);
end
%Plot MaxStress vs windspeed
Sigma_allowed=108; %Mpa
sigmafig=figure;
sigmaAx=axes(sigmafig); hold(sigmaAx,'on'); grid(sigmaAx,'on');
plot(sigmaAx,aUinf/alpha,MaxSigma*micro,'linewidth',2);
plot(sigmaAx,sigmaAx.XLim,[Sigma_allowed,Sigma_allowed],'r','linewidth',2);
plot(sigmaAx,[min(wr)/alpha,min(wr)/alpha],sigmaAx.YLim,'--m','linewidth',1)
plot(sigmaAx,[min(wr)/(3*alpha),min(wr)/(3*alpha)],sigmaAx.YLim,'--m',...
    'linewidth',1,'HandleVisibility','off')
plot(sigmaAx,[max(wr)/alpha,max(wr)/alpha],sigmaAx.YLim,'--g','linewidth',1)
plot(sigmaAx,[max(wr)/(3*alpha),max(wr)/(3*alpha)],sigmaAx.YLim,...
    '--g','linewidth',1,'HandleVisibility','off')
xlabel(sigmaAx,'U_\infty [m/s]'); 
ylabel(sigmaAx,'\sigma [Mpa]'); 
legend(sigmaAx,'Maximal stress in beam','Allowed stress',...
    '[0.33,1]*w_{r1}/\alpha','[0.33,1]*w_{r2}/\alpha');
%% B-initalize
M=9e3; %kg
d=3; %m
Theta=20*pi/180; %rad
cz=2e4; %Ns/m
%% B-rotate simulation
%Simulation for z kova
phi=@(v) atan2(v,l);
R=@(phi) [cos(phi),-sin(phi);sin(phi),cos(phi)];
zKova=[-sin(Theta);-cos(Theta)]; %[e1,e2]
RzKova=R(phi(0.3))*zKova;
AzKova=[sin(phi(0.3)-Theta);-cos(phi(0.3)-Theta)];
%Rotate in another 90 degrees to fit fked up axes chosen
zKova=R(pi/2)*zKova;
RzKova=R(pi/2)*RzKova;
AzKova=R(pi/2)*AzKova;
zKovafig=figure;
zKovaAx=axes(zKovafig);
hold(zKovaAx,'on'); grid(zKovaAx,'on'); axis(zKovaAx,'equal');
plot([0,zKova(1)],[0,zKova(2)]); 
plot([0,RzKova(1)],[0,RzKova(2)]);
plot([0,AzKova(1)],[0,AzKova(2)],'--');
%Simulation works!!!!!!!!
%% B-Q1alef
syms v u z dv du dz sEI sTheta sm sM sd sl sL skz sc scz real %d~time derivative
rm=[v,u];
drm=[dv,du];
rM=[v+(sd+z)*sin(atan2(v,sl)-sTheta),u-(sd+z)*cos(atan2(v,l)-sTheta)];
drM=[dv+(sd+z)*cos(atan2(v,sl)-sTheta)*(sl)/(v^2+sl^2)*dv+dz*sin(atan2(v,sl)-sTheta),...
    du+(sd+z)*sin(atan2(v,sl)-sTheta)*(sl)/(v^2+sl^2)*dv-dz*cos(atan2(v,sl)-sTheta)];
%% B-Q1bet
T=0.5*sm*(drm*drm')+0.5*sM*(drM*drM');
Mmat=simplify([diff(T,dv,2),diff(diff(T,dv),du),diff(diff(T,dv),dz);
        diff(diff(T,du),dv),diff(T,du,2),diff(diff(T,du),dz);
        diff(diff(T,dz),dv),diff(diff(T,dz),du),diff(T,dz,2)]);
u0=0; v0=0; z0=0;
M0=double(subs(Mmat,[u,v,z,sL,sl,sTheta,sM,sm,sd],[u0,v0,z0,L,l,Theta,M,m,d]));
%% B-Q1gimel
V=0.5*sEI*int(ddw^2,x,0,sL)+0.5*skz*z^2; %first V taken from A-Q1
Km=simplify([diff(V,v,2),diff(diff(V,v),u),diff(diff(V,v),z);
        diff(diff(V,u),v),diff(V,u,2),diff(diff(V,u),z);
        diff(diff(V,z),v),diff(diff(V,z),u),diff(V,z,2)]);
K0kz=simplify(subs(Km,[v,u,z,sL,sl,sEI],[v0,u0,z0,L,l,E*I])); %K0 - not complete doubled
%% B-Q1daled
D=0.5*sc*(drm*drm')+0.5*scz*(dz*dz');
Cmat=simplify([diff(D,dv,2),diff(diff(D,dv),du),diff(diff(D,dv),dz);
            diff(diff(D,du),dv),diff(D,du,2),diff(diff(D,du),dz);
            diff(diff(D,dz),dv),diff(D,dz,du),diff(D,dz,2)]);
du0=0; dv0=0; dz0=0;
C0=double(subs(Cmat,[u,v,z,du,dv,dz,sL,sl,sTheta,sM,sm,sd,sc,scz],[u0,v0,z0,du0,dv0,dz0,L,l,Theta,M,m,d,c,cz]));
%% B-Q2
%Taking a good guess for a working kz
kz=M*(wr_2(1,1)); %wr(1)~28rad/s
[MaxSigma,MaxZ]=deal(zeros(size(aUinf)));
K0=double(subs(K0kz,skz,kz));
K0=0.5*(K0+K0'); %enforce sym
[PHI,wr_3]=eig(K0,M0); %solves Av=lambda*Bv ~ K-(wn^2)*M
wr=sqrt(diag(wr_3)); %change wr to new system
Gama=PHI'*C0*PHI;
Lambda=PHI'*K0*PHI;
H=@(w) diag(1./diag(((-w^2)*eye(3)+(1i*w)*Gama+Lambda))); %3x3 matrix, normalized mass
for kk=1:length(aUinf) %loop around and build MaxSigma array
    eta1=(H(3*aUinf(kk))*PHI'*[B*exp(psi*1i);0;0]); %for P input
    eta2=(H(aUinf(kk))*PHI'*[0;A;0]); %for F input
    %Build time vector to sum upon
    Cycles=5;
    T=2*pi/(aUinf(kk));
    Fs=10*(1/T);
    t=0:1/Fs:(Cycles*T);
    %Calculate generalized co-ordinates q
    q=PHI*real((eta1*exp(3i*aUinf(kk)*t)+eta2*exp(1i*aUinf(kk)*t)));
    MomentMax=arrayfun(@(v,u) Momentx0(v,u),q(1,:),q(2,:));
    MaxSigma(kk)=max(sigma(MomentMax));
    MaxZ(kk)=max(abs(q(3,:)));
end

kzfig=figure;
%Plot max Z as a function of windspeed
MaxZAx=subplot(2,1,2,'parent',kzfig); hold(MaxZAx,'on'); grid(MaxZAx,'on');
ylabel(MaxZAx,'maximum |z| ordinate [m]'); 
xlabel(MaxZAx,'U_\infty [m/s]');
plot(MaxZAx,Uinf,MaxZ,'linew',2);

%Plot MaxStress vs windspeed
Sigma_allowed=108; %Mpa
sigmaAx=subplot(2,1,1,'parent',kzfig); hold(sigmaAx,'on'); grid(sigmaAx,'on');
plot(sigmaAx,Uinf,MaxSigma*micro,'linewidth',2);
plot(sigmaAx,sigmaAx.XLim,[Sigma_allowed,Sigma_allowed],'r','linewidth',2);
plot(sigmaAx,[wr(1)/alpha,wr(1)/alpha],sigmaAx.YLim,'--m','linewidth',1)
plot(sigmaAx,[wr(1)/(3*alpha),wr(1)/(3*alpha)],sigmaAx.YLim,'--m',...
    'linewidth',1,'HandleVisibility','off')
plot(sigmaAx,[wr(2)/alpha,wr(2)/alpha],sigmaAx.YLim,'--g','linewidth',1)
plot(sigmaAx,[wr(2)/(3*alpha),wr(2)/(3*alpha)],sigmaAx.YLim,...
    '--g','linewidth',1,'HandleVisibility','off')
plot(sigmaAx,[wr(3)/alpha,wr(3)/alpha],sigmaAx.YLim,'--b','linewidth',1)
plot(sigmaAx,[wr(3)/(3*alpha),wr(3)/(3*alpha)],sigmaAx.YLim,...
    '--b','linewidth',1,'HandleVisibility','off')
xlabel(sigmaAx,'U_\infty [m/s]'); 
ylabel(sigmaAx,'\sigma [Mpa]'); 
legend(sigmaAx,'Maximal stress in beam','Allowed stress',...
    '[1,3]*w_{r1}/\alpha','[1,3]*w_{r2}/\alpha','[1,3]*w_{r3}/\alpha');