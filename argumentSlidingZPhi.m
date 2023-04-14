%% Sliding Surface Model
function [sys,x0,str,ts] = argumentSlidingZPhi(t,x,u,flag) 
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
sizes.NumInputs = 12; 
sizes.DirFeedthrough = 4; 
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
hl  = u(12);

u0 = [u01 u02]';


Z11d   = 1 +t;

dZ12d  = 1;
ddZ12d = 0;
Phi11d   = 0.8;
dPhi12d  = 0;
ddPhi12d = 0;

Yd = [Z11d dZ12d Phi11d dPhi12d]'; %reference vector
dYd = [dZ12d ddZ12d dPhi12d ddPhi12d]'; %derivative of reference vector

%  ml   = 0.6;
%  hl   = 0.5;
% m0 = 0.6;
kc  = 0.08;  %mass flow
Sn  = 0.025;  %surface area of the nozzle
ax = 0.36;    %length of the tank
by = 0.36;    %breadth of the tank
rho  = 1002; %density of the liquid 
dm   = - kc*rho*Sn*hl;% rate of change of the mass
Mq = 12; %Quadrotor Mass = 12kg
M    = Mq + ml; cz = 0.03; %cz is the distance b/n the center of masses 


C0 = ml/12;

Ilx = by.^2 + hl.^2;


Il = C0.*Ilx; %Inertia of the liquid

Iq = diag([1.5 1.5 3]);
Ix = Iq(1,1) + Il;

PI1 = -(dm*cz^2)/Ix - dm*( by^2 + 3*ml^2)/(12*ax^2 * by^2 * rho^2 *Ix); 

 
f = [-(dm/M)*x12; PI1*x22];
 
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
     0 1/Ix];

 
uo = u0;
     
Iz = Hz*(A*z + Fz + B*g0*uo);

sys(1) = Iz(1);
sys(2) = Iz(2);
sys(3) = Iz(3);
sys(4) = Iz(4);





