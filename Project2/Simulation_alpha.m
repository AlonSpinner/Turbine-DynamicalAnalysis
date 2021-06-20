%% Simulation alpha
%choose parametrs
aUinf=28; %rad/s
kz=1e5; %N/m
vfactor=1;
%% initial parameters
M=9e3; %kg
d=3; %m
Theta=20*pi/180; %rad
mili=1e-3; micro=1e-6;
E=210e9; %Pa
L=65; %m
l=2; %m
m=90e3; %kg
A=1.9e6; %N
B=0.8e6; %N;
alpha=10; %rad/m
psi=pi/3; %rad

%drawing parameters
ballr=0.5;
tc=linspace(0,2*pi,20);
SpringR=0.2; SpringN=1;
spr=Spring(SpringR,SpringN);
%% calculate q
%Basic matrices and modal system
M0=1e5*[1.446217007612195,  -0.046172719348965,  -0.030781812899310;
  -0.046172719348965,   0.990000000000000,  -0.084572335870732;
  -0.030781812899310,  -0.084572335870732,   0.090000000000000];
C0=1e5*[3.413815572471545  -0.102606042997701  -0.068404028665134;
  -0.102606042997701   2.400000000000000  -0.187938524157182;
  -0.068404028665134  -0.187938524157182   0.200000000000000];
K0=[7082909133077403/66560,21248727399232209/4326400,0;
    21248727399232209/4326400,21248727399232209/70304000,0;
    0,0,kz];
K0=0.5*(K0+K0'); %enforce sym
[PHI,wr_3]=eig(K0,M0); %solves Av=lambda*Bv ~ K-(wn^2)*M
wr=sqrt(diag(wr_3)); %change wr to new system
Gama=PHI'*C0*PHI;
Lambda=PHI'*K0*PHI;

%time response
H=@(w) diag(1./diag(((-w^2)*eye(3)+(1i*w)*Gama+Lambda))); %3x3 matrix, normalized mass
eta1=(H(3*aUinf)*PHI'*[B*exp(psi*1i);0;0]); %for P input
eta2=(H(aUinf)*PHI'*[0;A;0]); %for F input
%Build time vector to sum upon
Cycles=20;
T=min(2*pi/(aUinf),10);
Fs=10*(1/T);
t=0:1/Fs:(Cycles*T);
%Calculate generalized co-ordinates q
q=PHI*real((eta1*exp(3i*aUinf*t)+eta2*exp(1i*aUinf*t)));
%% simulation
SimFig=figure;
SimAx=axes(SimFig);
hold(SimAx,'on'); grid(SimAx,'on');
axis(SimAx,'equal'); %SimAx.XLimMode='manual';
xlim(SimAx,[-3,7]); ylim(SimAx,[-3,3]);
xlabel(SimAx,'[m]'); ylabel(SimAx,'[m]');
rM=@(q) [q(1)+(d+q(3))*sin(atan2(q(1),l)-Theta),q(2)-(d+q(3))*cos(atan2(q(1),l)-Theta)];
rm=@(q) [q(1),q(2)];
w=@(q) (4225*(3*q(2)+L*atan2(q(1),l)))/L^2-(274625*(2*q(2)+L*atan2(q(1),l)))/L^3;
R=@(phi) [cos(phi),-sin(phi);sin(phi),cos(phi)];
l0=[2;0];
d0=3*[cos(Theta);-sin(Theta)];
Colors=lines(4);
for kk=1:length(t)
cla(SimAx);
qloop=q;
qloop(1,:)=qloop(1,:)*vfactor;
%Obtain CG points and rotate them to fit matlab axes
cgm=rm(qloop(:,kk));
cgM=rM(qloop(:,kk));
cgM=R(pi/2)*cgM';
cgm=R(pi/2)*cgm';
%calc l and d links in time
lfv=R(atan2(qloop(1,kk),l))*l0;
dfv=R(atan2(qloop(1,kk),l))*d0;
%draw link l
plot(SimAx,cgm(1)+[0,lfv(1)],cgm(2)+[0,lfv(2)],'color',Colors(3,:),'linew',5);
%Draw masses
patch(SimAx,cgm(1)+ballr*cos(tc),cgm(2)+ballr*sin(tc),Colors(1,:));
patch(SimAx,cgM(1)+0.5*ballr*cos(tc),cgM(2)+0.5*ballr*sin(tc),Colors(2,:));
%draw link d
plot(SimAx,cgm(1)+[0,dfv(1)],cgm(2)+[0,dfv(2)],'color',Colors(4,:),'linew',5);
%draw kz spring
[x,y]=spr.getSpr(cgM,dfv+cgm);
plot(SimAx,x,y,'k','linew',2);

pause(0.1)
end