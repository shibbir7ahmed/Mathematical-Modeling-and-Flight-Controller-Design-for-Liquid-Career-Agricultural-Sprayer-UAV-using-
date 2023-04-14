%% Sliding Surface Model
function [sys,x0,str,ts] = argumentSlidingXY(t,x,u,flag) 
switch flag,
    
    case 0, [sys,x0,str,ts] = mdlInitializeSizes;
    
    case 1, sys = mdlDerivatives(t,x,u);
    
    case 3, sys = mdlOutputs(t,x,u);
    
    case {2, 4, 9 } 
        
        sys = [];

    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts] = mdlInitializeSizes
sizes = simsizes;
sizes.NumContStates = 0;
sizes.NumDiscStates = 0;
sizes.NumOutputs = 4; 
sizes.NumInputs = 11; 
sizes.DirFeedthrough = 2; 
sizes.NumSampleTimes = 0; 
sys = simsizes(sizes); 
x0 = [];
str = [];
ts = [];


function sys=mdlOutputs(t,x,u)
u01 = u(1);
u02 = u(2);
z1  = u(3);
z2  = u(4);
z3  = u(5);
z4  = u(6);
x11 = u(7);
x12 = u(8);
x21 = u(9);
x22 = u(10);
ml  = u(11);

u0 = [u01 u02]';


X11d = cos(0.5*t);

if t >= 0 && t < 10
    
        Y11d = 3; %sin(t) + sin(0.5*t);
end

if t >= 10 && t < 20
    
        Y11d = 1; %sin(t) + sin(0.5*t);
end

if t >= 20 
    
        Y11d = 3; %sin(t) + sin(0.5*t);
end

dX12d  = -0.5*sin(0.5*t);
dY12d  = 0;%cos(t) + 0.5*cos(0.5*t);

ddX12d = -(0.5)*(0.5)*sin(0.5*t);
ddY12d = 0;

Yd = [X11d dX12d Y11d dY12d]'; %reference vector


dYd = [dX12d ddX12d dY12d ddY12d]'; %derivative of reference vector

% ml = 0.6;
Mq = 12; %Quadrotor Mass = 12kg
dml = 0.8;
M = ml + Mq; %J2 = 1; m1 =2; m2 = 2; r = 0.1;l = 0.5; b = 0.4; g = 9.8;
 
f = [-(dml/M)*x12; -(dml/M)*x22];
 
z  = [z1 z2 z3 z4]';

 A = [0 1 0 0;
     0 0 0 0;        
     0 0 0 1;        
     0 0 0 0];

B = [0 0;
     1 0;
     0 0;
     0 1];

 Fz = B*f - (dYd - A*Yd);
 
Hz = eye(4);
 

g0 = [1/M 0;
     0 1/M];

 
uo = u0;
     
Iz = Hz*(A*z + Fz + B*g0*uo);

sys(1) = Iz(1);
sys(2) = Iz(2);
sys(3) = Iz(3);
sys(4) = Iz(4);





