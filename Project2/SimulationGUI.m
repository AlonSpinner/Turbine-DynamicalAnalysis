%Vibration Project2 Result Simulation GUI 
%by Alon Spinner 305184335
function varargout = SimulationGUI(varargin)
% SimulationGUI MATLAB code for SimulationGUI.fig
%      SimulationGUI, by itself, creates a new SimulationGUI or raises the existing
%      singleton*.
%
%      H = SimulationGUI returns the handle to a new SimulationGUI or the handle to
%      the existing singleton*.
%
%      SimulationGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SimulationGUI.M with the given input arguments.
%
%      SimulationGUI('Property','Value',...) creates a new SimulationGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SimulationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SimulationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SimulationGUI

% Last Modified by GUIDE v2.5 06-Jan-2019 10:54:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SimulationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SimulationGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function SimulationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SimulationGUI (see VARARGIN)

% Choose default command line output for SimulationGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
function varargout = SimulationGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%% CallBacks
function aUinfEdit_Callback(hObject, eventdata, handles)
StopSimh=handles.StopSimToggle;
StopSimh.Value=1; %Stops any current simulations
function kzEdit_Callback(hObject, eventdata, handles)
StopSimh=handles.StopSimToggle;
StopSimh.Value=1; %Stops any current simulations
function vampEdit_Callback(hObject, eventdata, handles)
StopSimh=handles.StopSimToggle;
StopSimh.Value=1; %Stops any current simulations
function DamperToggle_Callback(hObject, eventdata, handles)
StopSimh=handles.StopSimToggle;
StopSimh.Value=1; %Stops any current simulations
function PauseEdit_Callback(hObject, eventdata, handles)
StopSimh=handles.StopSimToggle;
StopSimh.Value=1; %Stops any current simulations

function RunPush_Callback(hObject, eventdata, handles)
%Obtain user input and check it
aUinf=str2num(handles.aUinfEdit.String);
kz=str2num(handles.kzEdit.String);
vamp=str2num(handles.vampEdit.String);
PauseTime=handles.PauseEdit.String; %real or a number in string format
SimAx=handles.SimAx;
StopSimh=handles.StopSimToggle;
StopSimh.Value=1; %Stops any current simulations
StopSimh.Value=0; %allows to continue with new simulation

DamperToggle=handles.DamperToggle.Value;
if isempty(aUinf)||isempty(kz)||isempty(vamp)
    errordlg('Please check you inputs','A-lon');
    return
end

%ready SimAx
cla(SimAx,'reset');
hold(SimAx,'on'); grid(SimAx,'on');
axis(SimAx,'equal'); %SimAx.XLimMode='manual';
xlim(SimAx,[-3,4]); ylim(SimAx,[-3,3]);
xlabel(SimAx,'[m]'); ylabel(SimAx,'[m]');

%run simulation
if DamperToggle
    SimWMassDamper(StopSimh,SimAx,aUinf,kz,vamp,PauseTime)
else
    Sim(StopSimh,SimAx,aUinf,vamp,PauseTime)
end
function UpPush_Callback(hObject, eventdata, handles)
%Add +1 to the number of lines
vamp=str2num(handles.vampEdit.String);
if isempty(vamp)
    errordlg('Please check you inputs','A-lon');
    return
end
vamp=vamp+1;
handles.vampEdit.String=num2str(vamp);

StopSimh=handles.StopSimToggle;
StopSimh.Value=1; %Stops any current simulations
StopSimh.Value=0; %allows to continue with new simulation
function DownPush_Callback(hObject, eventdata, handles)
vamp=str2num(handles.vampEdit.String);
if isempty(vamp)
    errordlg('Please check you inputs','A-lon');
    return
end
vamp=vamp-1;
handles.vampEdit.String=num2str(vamp);

StopSimh=handles.StopSimToggle;
StopSimh.Value=1; %Stops any current simulations
StopSimh.Value=0; %allows to continue with new simulation
%% functions
function SimWMassDamper(StopSimh,SimAx,aUinf,kz,vfactor,PauseTime)
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
r_out=3; %m
r_in=2.5; %m
I=pi*(r_out^4-r_in^4)/(4); %m^4

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
C0= [220000,0,0;
    0,220000,0;
    0,0,20000];
K0=[7082909133077403/66560,21248727399232209/4326400,0;
    21248727399232209/4326400,21248727399232209/70304000,0;
    0,0,kz];
[PHI,~]=eig(K0,M0); %solves Av=lambda*Bv ~ K-(wn^2)*M
Gama=PHI'*C0*PHI;
Lambda=PHI'*K0*PHI;

%time response
H=@(w) diag(1./diag(((-w^2)*eye(3)+(1i*w)*Gama+Lambda))); %3x3 matrix, normalized mass
eta1=(H(3*aUinf)*PHI'*[B*exp(psi*1i);0;0]); %for P input
eta2=(H(aUinf)*PHI'*[0;A;0]); %for F input
%Build time vector to sum upon and play upon
Cycles=20;
T=min(2*pi/(aUinf),10);
Fs=20*(1/T);
t=0:(1/Fs):(Cycles*T);
%Calculate generalized co-ordinates q
q=PHI*real((eta1*exp(3i*aUinf*t)+eta2*exp(1i*aUinf*t)));
qloop=q;
qloop(1,:)=qloop(1,:)*vfactor;
%% simulation
rM=@(q) [q(1)+(d+q(3))*sin(atan2(q(1),l)-Theta),q(2)-(d+q(3))*cos(atan2(q(1),l)-Theta)];
rm=@(q) [q(1),q(2)];
R=@(phi) [cos(phi),-sin(phi);sin(phi),cos(phi)];
l0=[2;0];
d0=3*[cos(Theta);-sin(Theta)];
Colors=lines(4);
if strcmpi(PauseTime,'Real')
    PauseTime=(1/Fs);
else
    PauseTime=str2num(PauseTime);
end
for kk=1:length(t)
    if StopSimh.Value, return, end
    cla(SimAx);
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
    patch(SimAx,cgm(1)+ballr*cos(tc),cgm(2)+ballr*sin(tc),Colors(1,:),'edgec','n');
    patch(SimAx,cgM(1)+0.5*ballr*cos(tc),cgM(2)+0.5*ballr*sin(tc),Colors(2,:),'edgec','n');
    %draw link d
    plot(SimAx,cgm(1)+[0,dfv(1)],cgm(2)+[0,dfv(2)],'color',Colors(4,:),'linew',5);
    %draw kz spring
    [x,y]=spr.getSpr(cgM,dfv+cgm);
    plot(SimAx,x,y,'k','linew',2);
    
    pause(PauseTime)
end
function Sim(StopSimh,SimAx,aUinf,vfactor,PauseTime)
%% initial parameters
mili=1e-3; micro=1e-6;
E=210e9; %Pa
L=65; %m
l=2; %m
m=90e3; %kg
A=1.9e6; %N
B=0.8e6; %N;
alpha=10; %rad/m
psi=pi/3; %rad
r_out=3; %m
r_in=2.5; %m
I=pi*(r_out^4-r_in^4)/(4); %m^4

%drawing parameters
ballr=0.5;
tc=linspace(0,2*pi,20);
%% calculate q
%Basic matrices and modal system
Gama=[2.444444444444444,0;
      0   2.444444444444444];
Lambda=1e6*[0.000837771522979,0;
        0,1.184897116204642];
PHI=[0.000153791292281,-0.003329783679089;
  -0.003329783679089,-0.000153791292281];

%time response
H=@(w) diag(1./diag((-w^2)*eye(2)+(1i*w)*Gama+Lambda)); %2x2 matrix, %normalized mass
eta1=(H(3*aUinf)*PHI'*[B*exp(psi*1i);0]); %for P input
eta2=(H(aUinf)*PHI'*[0;A]); %for F input
%Build time vector to sum upon and play upon
Cycles=20;
T=min(2*pi/(aUinf),10);
Fs=20*(1/T);
t=0:(1/Fs):(Cycles*T);
%Calculate generalized co-ordinates q
q=PHI*real((eta1*exp(3i*aUinf*t)+eta2*exp(1i*aUinf*t)));
qloop=q;
qloop(1,:)=qloop(1,:)*vfactor;
%% simulation
rm=@(q) [q(1),q(2)];
R=@(phi) [cos(phi),-sin(phi);sin(phi),cos(phi)];
l0=[2;0];
Colors=lines(4);
if strcmpi(PauseTime,'Real')
    PauseTime=(1/Fs);
else
    PauseTime=str2num(PauseTime);
end
for kk=1:length(t)
    if StopSimh.Value, return, end
    cla(SimAx);
    %Obtain CG points and rotate them to fit matlab axes
    cgm=rm(qloop(:,kk));
    cgm=R(pi/2)*cgm';
    %calc l link in time
    lfv=R(atan2(qloop(1,kk),l))*l0;
    %draw link l
    plot(SimAx,cgm(1)+[0,lfv(1)],cgm(2)+[0,lfv(2)],'color',Colors(3,:),'linew',5);
    %Draw mass
    patch(SimAx,cgm(1)+ballr*cos(tc),cgm(2)+ballr*sin(tc),Colors(1,:),'edgec','n');
    
    pause(PauseTime)
end
%% Create Fcns
function aUinfEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aUinfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kzEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kzEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function vampEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function PauseEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function StopSimToggle_Callback(hObject, eventdata, handles)