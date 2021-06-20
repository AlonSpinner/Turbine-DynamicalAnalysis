%for n>3 eig returns negative frequencies and non-mass-normalized modal
%vectors

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

%choose n <----------------------------------------------------
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
K=blkdiag(Km,Kb); %is symmetrical for every n (checked)

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
M=blkdiag(Mm+Mbm,Mb)+Mpm; %is symmetrical for every n (checked)

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

%ensure symmatrey
Kc=0.5*(Kc+Kc');
Mc=0.5*(Mc+Mc');

[PHI,wr2]=eig(Kc,Mc);

disp('Squared Natrual Frequencies');
disp(diag(wr2)');
disp('Modal Masses');
disp(diag(PHI'*Mc*PHI)');