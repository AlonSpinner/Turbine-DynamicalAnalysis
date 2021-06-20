%% Initalize 
syms  xi chi real
L=65; %m
l=10; %m
rb=1; %m
pm=90e3;  %kg %mass of point mass
rho=8e3; %kg/m^3
E=210e9; %Pa
rm_in=2.5-xi/L;  %m, symbolic
rm_out=3-xi/L; %m, symbolic

%Calculating area and moment of inertia
Am=pi*(rm_out^2 -rm_in^2); %m^2
Ab=pi*rb^2; %m^2
Im=pi/4*(rm_out^4-rm_in^4); %m^2
Ib=pi/4*rb^4; %m^2

%% Q4 Calculating matrices 
% w=Psy'*q; %mast deflection
% u=Phi'*p %beam deflection
% z=[q;p]; - general ordinates

%Number of base functions for rayligh ritz
n=2;

Psy=Chebypoly(xi,n,0,L);
dPsy=diff(Psy,xi);
ddPsy=diff(dPsy,xi);

Phi=Chebypoly(chi,n,0,l);
dPhi=diff(Phi,chi);
ddPhi=diff(dPhi,chi);

%Stiffness matrix
Km=double(E*int(Im*ddPsy*ddPsy',0,L)); %mast stiffness matrix
Kb=double(E*int(Ib*ddPhi*ddPhi',0,l)); %beam stiffness matrix
K=blkdiag(Km,Kb);

%Mass matrix
PsyL=subs(Psy,xi,L); %prep up for calculations
dPsyL=subs(dPsy,xi,L);
Phil=subs(Phi,chi,l);
bm=rho*Ab*l; %mass of beam
Mm=double(rho*int(Am*Psy*Psy',0,L)); %mast mass matrix
Mbm=double(bm*(PsyL*PsyL'+(l^2)/3*dPsyL*dPsyL')); %mast mass matrix - beam cg motion
Mb=double(rho*int(Ab*Phi*Phi',0,l)); %beam mass matrix
Mpm=pm*double([PsyL*PsyL'+(l^2)*dPsyL*dPsyL', l*dPsyL*Phil'; %points mass matrix
   l*Phil*dPsyL',Phil*Phil']);
M=blkdiag(Mm+Mbm,Mb)+Mpm;

%Defining constraints - cantilever fixed end
c = double([
            [subs(Psy',xi,0),0*Phi']; %w(0)=0
            [subs(dPsy',xi,0),0*dPhi'];%w'(0)=0
            [0*Psy',subs(Phi',chi,0)];%u(0)=0
            [0*dPsy',subs(dPhi',chi,0)];%u'(0)=0
            ]);
C=null(c,'r');

%Apply constraints on matrices
Kc=C'*K*C;
Mc=C'*M*C;
%% Q5 calcualting generalized forces

syms Uinf t F A alpha real 
% A=35*10^4; %N/m
% alpha=3.6*3.5; %(rad/m)
F=A*cos(alpha*Uinf*t);
Q=int(F*[Psy;0*Phi],0,L);
Qc=C'*Q;

%% Q6 selection of number of modes 
nmax=20; %maximum n - number of polynomials for each deflection function
wrMat=nan(2*(nmax+1)-4,nmax); %Initalize.
for n=2:nmax
    
Psy=Chebypoly(xi,n,0,L);
dPsy=diff(Psy,xi);
ddPsy=diff(dPsy,xi);

Phi=Chebypoly(chi,n,0,l);
dPhi=diff(Phi,chi);
ddPhi=diff(dPhi,chi);

%Stiffness matrix
Km=double(E*int(Im*ddPsy*ddPsy',0,L)); %mast stiffness matrix
Kb=double(E*int(Ib*ddPhi*ddPhi',0,l)); %beam stiffness matrix
K=blkdiag(Km,Kb);

%Mass matrix
PsyL=subs(Psy,xi,L); %prep up for calculations
dPsyL=subs(dPsy,xi,L);
Phil=subs(Phi,chi,l);
bm=rho*Ab*l; %mass of beam
Mm=double(rho*int(Am*Psy*Psy',0,L)); %mast mass matrix
Mbm=double(bm*(PsyL*PsyL'+(l^2)/3*dPsyL*dPsyL')); %mast mass matrix - beam cg motion
Mb=double(rho*int(Ab*Phi*Phi',0,l)); %beam mass matrix
Mpm=pm*double([PsyL*PsyL'+(l^2)*dPsyL*dPsyL', l*dPsyL*Phil'; %points mass matrix
   l*Phil*dPsyL',Phil*Phil']);
M=blkdiag(Mm+Mbm,Mb)+Mpm;

%Defining constraints - cantilever fixed end
c = double([
            [subs(Psy',xi,0),0*Phi']; %w(0)=0
            [subs(dPsy',xi,0),0*dPhi'];%w'(0)=0
            [0*Psy',subs(Phi',chi,0)];%u(0)=0
            [0*dPsy',subs(dPhi',chi,0)];%u'(0)=0
            ]);
C=null(c);

%Apply constraints on matrices
Kc=C'*K*C;
Mc=C'*M*C;

%ensure symmatrey
Kc=0.5*(Kc+Kc');
Mc=0.5*(Mc+Mc');

[PHI,wr2]=eig(Kc,Mc); %solves Av=lambda*Bv ~ K-(wn^2)*M
wr=sqrt(diag(wr2));
wrMat(1:length(wr),n)=wr;
end

%decide on number of modes of importance
N=6;

%plot natrual frequencies
[SwrMat,Indx]=sort(wrMat); %sort 
TSwrMat=SwrMat'; %transpose
WrFig=figure('Name','WrFig');
WrAx=axes(WrFig);
hold(WrAx,'on'); grid(WrAx,'on');
plot(WrAx,1:nmax,TSwrMat(:,1:N),'--o','MarkerSize',5);
xlabel(WrAx,'n - base polynomial number');
ylabel(WrAx,'Natrual Frequency [rad/s]');
lgndWn=repmat({'w_n '},N,1);
lgndNums=arrayfun(@num2str,[1:N],'un',0);
lgnd=strcat(lgndWn,lgndNums');
legend(WrAx,lgnd);
%% plot modes
Modes=C*PHI;
SModes=Modes(:,Indx(:,end)); %sorted
for k=1:N
    Plot_Turbine(5000*SModes(:,k),k+1); %figure number k is reserved for WrFig
end

%% Simulation for Q6
%Mode Simulation
k=1;
SimFig=figure;
SimAx=axes(SimFig);
SimMode=5000*SModes(:,k);
wn=TSwrMat(end,k); %rad/s
T=2*pi/wn;
t=linspace(0,2*T,100);
SimMode_t=real(SimMode*exp(1i*wn*t));
for i=1:length(t)
Plot_Turbine(SimMode_t(:,i),1); %figure number k is reserved for WrFig
drawnow
pause(0.1); %seconds
end

%% Q7
%Solve with n=8, but take only first 5 modes for solution
n=8;

Psy=Chebypoly(xi,n,0,L);
dPsy=diff(Psy,xi);
ddPsy=diff(dPsy,xi);

Phi=Chebypoly(chi,n,0,l);
dPhi=diff(Phi,chi);
ddPhi=diff(dPhi,chi);

%Stiffness matrix
Km=double(E*int(Im*ddPsy*ddPsy',0,L)); %mast stiffness matrix
Kb=double(E*int(Ib*ddPhi*ddPhi',0,l)); %beam stiffness matrix
K=blkdiag(Km,Kb);

%Mass matrix
PsyL=subs(Psy,xi,L); %prep up for calculations
dPsyL=subs(dPsy,xi,L);
Phil=subs(Phi,chi,l);
bm=rho*Ab*l; %mass of beam
Mm=double(rho*int(Am*Psy*Psy',0,L)); %mast mass matrix
Mbm=double(bm*(PsyL*PsyL'+(l^2)/3*dPsyL*dPsyL')); %mast mass matrix - beam cg motion
Mb=double(rho*int(Ab*Phi*Phi',0,l)); %beam mass matrix
Mpm=pm*double([PsyL*PsyL'+(l^2)*dPsyL*dPsyL', l*dPsyL*Phil'; %points mass matrix
   l*Phil*dPsyL',Phil*Phil']);
M=blkdiag(Mm+Mbm,Mb)+Mpm;

%Defining constraints - cantilever fixed end
c = double([
            [subs(Psy',xi,0),0*Phi']; %w(0)=0
            [subs(dPsy',xi,0),0*dPhi'];%w'(0)=0
            [0*Psy',subs(Phi',chi,0)];%u(0)=0
            [0*dPsy',subs(dPhi',chi,0)];%u'(0)=0
            ]);
C=null(c);

%Apply constraints on matrices
Kc=C'*K*C;
Mc=C'*M*C;

%ensure symmatrey
Kc=0.5*(Kc+Kc');
Mc=0.5*(Mc+Mc');

[PHI,wr2]=eig(Kc,Mc);
wr=sqrt(diag(wr2));
[Swr,Indx]=sort(wr); %sort
Swr=abs(Swr'); %abs or real?
SPHI=PHI(:,Indx); %PHI to work with

zeta=0.03;
Gama=diag(2*zeta*Swr);
Lambda=diag(Swr.^2);

A=35e4; %N/m
alpha=3.5; %(rad*H/(km*s))
F0=A;
Q0=double(int(F0*[Psy;0*Phi],xi,0,L));

%Initalize vectors
ResUinf=500;
Uinf=linspace(0,60,ResUinf);
aUinf=alpha*Uinf;
ResXiChi=20;
xiVec=linspace(0,L,ResXiChi);
chiVec=linspace(0,l,ResXiChi);
[SmaxMast,SmaxBeam,xiMax,chiMax]=deal(zeros(size(aUinf)));

%initalize functions
H=@(w) diag(1./diag((-w^2)*eye(2*(n-1))+(1i*w)*Gama+Lambda)); %2*(n-1) is length of z
SigmaMast=@(ddw,xiVal) double(subs(rm_out*E*ddw,xi,xiVal)); %Pa
SigmaBeam=@(ddu,chiVal) double(subs(rb*E*ddu,chi,chiVal)); %Pa

%number of modes to calculate by
N=5;

%initalize time
Cycles=2;
T=2*pi/(Swr(1)); %first natural frequency
Fs=10*(Swr(N))/(2*pi);
t=0:1/Fs:(Cycles*T);
[SmaxMast_tt,SmaxBeam_tt,xiInd_tt,chiInd_tt]=deal(zeros(size(t)));

tic;
for kk=1:length(aUinf) %loop around and build MaxSigma array
    Htf=H(aUinf(kk));
    Htf(N+1:end,N+1:end)=0; %zero the transfer function for non relevant modes
    z0=C*SPHI*Htf*SPHI'*C'*Q0;
    %Build time vector to sum upon
    z=real(z0*exp(1i*aUinf(kk)*t)); %Calculate generalized co-ordinates q in time 
    q=z(1:n+1,:); %(n+1)Xlength(t) matrix
    p=z(n+2:end,:); %(n+1)Xlength(t) matrix
    for tt=1:length(t)
        ddw=ddPsy'*q(:,tt);
        ddu=ddPhi'*p(:,tt);
        SMast=SigmaMast(ddw,xiVec);
        SBeam=SigmaBeam(ddu,chiVec);
        [SmaxMast_tt(tt),xiInd_tt(tt)]=max(abs(SMast));
        [SmaxBeam_tt(tt),chiInd_tt(tt)]=max(abs(SBeam));
    end
    [SmaxMast(kk),MastInd]=max(SmaxMast_tt);
    [SmaxBeam(kk),BeamInd]=max(SmaxBeam_tt);
    xiMax(kk)=xiVec(xiInd_tt(MastInd));
    chiMax(kk)=chiVec(chiInd_tt(BeamInd));
    disp([num2str(kk),'/',num2str(length(aUinf))]);
    toc;
end
%Plot MaxStress vs windspeed
%Plot MaxStress vs windspeed
micro=1e-6;
Sigma_allowed=108; %Mpa
sigmafig=figure;
sigmaAx=axes(sigmafig); hold(sigmaAx,'on'); grid(sigmaAx,'on');
xlabel(sigmaAx,'U_\infty [km/h]'); 
ylabel(sigmaAx,'\sigma [Mpa]'); 
plot(sigmaAx,Uinf,SmaxMast*micro,'linewidth',2);
plot(sigmaAx,Uinf,SmaxBeam*micro,'linewidth',2);
plot(sigmaAx,xlim(sigmaAx),Sigma_allowed*[1,1],'linewidth',2);
legend(sigmaAx,'Maximal stress in Mast','Maximal stress in Beam','Allowed stress')
%Plot xiMax and chiMax vs windspeed
Xfig=figure;
XAx=axes(Xfig); hold(XAx,'on'); grid(XAx,'on');
xlabel(XAx,'U_\infty [km/h]'); 
ylabel(XAx,'\xi/\chi [m]');
plot(XAx,Uinf,xiMax,'linewidth',2);
plot(XAx,Uinf,chiMax,'linewidth',2);
legend(XAx,'xi of maximum stress in mast','chi of maximum stress in beam')