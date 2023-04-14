%% Sliding Surface Model
function [sys,x0,str,ts] = argumentSlidingThetaPsi(t,x,u,flag) 
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


Theta11d   = sin(t)+sin(0.5*t);
dTheta12d  = cos(t)+(0.5)*cos(0.5*t);
ddTheta12d = -sin(t)-(0.5)*(0.5)*sin(0.5*t);

% Y11d   = cos(t)+cos(t);
% dY12d  = -sin(t)-sin(t);
% ddY12d = -cos(t)-cos(t);

Psi11d   = sin(t)+sin(0.2*t);
dPsi12d  = cos(t)+(0.2)*cos(0.2*t);
ddPsi12d = -sin(t)-(0.2)*(0.2)*sin(0.5*t);

Yd = [Theta11d dTheta12d Psi11d dPsi12d]'; %reference vector
dYd = [dTheta12d ddTheta12d dPsi12d ddPsi12d]'; %derivative of reference vector

% ml   = 0.6;
% hl   = 0.5;
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
Ily = ax.^2+hl.^2;
Ilz = ax.^2+by.^2;

Il = C0.*[Ilx      0           0;
           0      Ily          0;
           0       0          Ilz]; %Inertia of the liquid

Iq = diag([1.5 1.5 3]); 
Iy = Iq(2,2) + Il(2,2);
Iz = Iq(3,3) + Il(3,3);


PI2 = -(dm*cz^2)/Iy - dm*( ax^2 + 3*ml^2)/(12*ax^2 * by^2 * rho^2 *Iy);
PI3 = - dm*( ax^2 + by^2)/(12*Iz); 
 
f = [PI2*x12; PI3*x22];
 
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
 

g0 = [1/Iy    0;
      0      1/Iz];

 
uo = u0;
     
Iz = Hz*(A*z + Fz + B*g0*uo);

sys(1) = Iz(1);
sys(2) = Iz(2);
sys(3) = Iz(3);
sys(4) = Iz(4);





