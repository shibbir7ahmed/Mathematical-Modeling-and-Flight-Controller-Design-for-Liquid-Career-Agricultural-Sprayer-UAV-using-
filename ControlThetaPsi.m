function [sys,x0,str,ts] = ControlThetaPsi(t,x,u,flag) 
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
global Ph
sizes = simsizes;
sizes.NumContStates = 75;
sizes.NumDiscStates = 0;
sizes.NumOutputs = 16; 
sizes.NumInputs = 16; 
sizes.DirFeedthrough = 4; 
sizes.NumSampleTimes = 0; 
sys = simsizes(sizes); 


y1d0   = sin(0)+sin(0.5*(0));
dy1d0  = cos(0)+(0.5)*cos(0.5*(0));
y2d0   = sin(0)+sin(0.2*(0));
dy2d0  = cos(0)+(0.2)*cos(0.2*(0));

X0  =[0.2, -0.2, 2, 2]';




Yd0 = [y1d0 dy1d0 y2d0 dy2d0]';
Eta = X0 -Yd0;


R = .10;%initializing the disturbance compensation gain of Eq. 20

% Wf=-1+(1+1)*rand(11,1);
% Wf=   [-0.5869
%    -0.2622
%     0.9410
%    -0.9271
%     0.5210
%    -0.8052
%    -0.5223
%     0.5088
%    -0.5685
%    -0.6398
%     0.7947]';

Wf = [0.5829
    0.8996
    0.5070
   -0.9630
   -0.3803
   -0.0061
   -0.3639
   -0.4310
   -0.8662
    0.1038
   -0.2719]';

Ws =Wf;
% Ws  = [0.0959
%     0.3251
%     0.8846
%    -0.6310
%    -0.5141
%    -0.9973
%     0.9904
%    -0.5750
%     0.9052
%    -0.0149
%     0.6978]';%second weight for RBF (11 by 1) of Eq.23


Ph  = [0 0 0 0 0 0 0 0 0 0 0]';%Initializing RBF (11 by 1) of Eq. 23


V0 = -1 + (1+1)*rand(11,4);


Theta = 0.1*[0 0;0 0];  %===========================================
           

for i =1:1:4
    x0(i) = Eta(i);
end

for ii =1:1:1
    x0(4+ii) = R(ii);
end

for iii = 1:1:11
    x0(5+iii) = Wf(iii);
end

for ij = 1:1:11
    x0(16+ij) = Ws(ij);
end

for i = 1:1:2
    x0(27+i) = Theta(1,i);
end

for i = 1:1:2
    x0(29+i) = Theta(2,i);
end

for i = 1:1:11
    x0(31+i) = V0(i,1);
end

for i = 1:1:11
    x0(42+i) = V0(i,2);
end

for i = 1:1:11
    x0(53+i) = V0(i,3);
end

for i = 1:1:11
    x0(64+i) = V0(i,4);
end


str = [];
ts = [];

function sys=mdlDerivatives(t,x,u)
global Ph

z1  = u(1);
z2  = u(2);
z3  = u(3);
z4  = u(4);
x11 = u(5);
x12 = u(6);
x21 = u(7);
x22 = u(8);
S1  = u(9);
S2  = u(10);
S3  = u(11);
S4  = u(12);
Et1 = u(13);
Et2 = u(14);
ml  = u(15);
hl  = u(16);
% Dz1 = u(15);
% Dz2 = u(16);
% Dz3 = u(17);
% Dz4 = u(18);
% 
% Dz = [Dz1 Dz2 Dz3 Dz4]';

Et = [Et1 Et2]';

S = [S1 S2 S3 S4]';

X = [x11 x12 x21 x22]';


for li =1:1:4
    Eta(li)=x(li);
end

for il =1:1:1
    R(il) = x(4+il);
end

for ir = 1:1:11  
    Wf(ir) = x(5+ir);
end

for ri = 1:1:11
    Ws(ri) = x(16+ri);
end

W = [Wf' Ws'];

for i = 1:1:2
    Theta(1,i) = x(27+i);
end

for i = 1:1:2
    Theta(2,i) = x(29+i);
end

for i = 1:1:11
    V0(i,1) = x(31+i);
end

for i = 1:1:11
    V0(i,2) = x(42+i);
end

for i = 1:1:11
    V0(i,3) = x(53+i);
end

for i = 1:1:11
    V0(i,4) = x(64+i);
end


Z = [z1 z2 z3 z4]';

b = 0.1;

C = 2*[-2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2;
     -2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2;
     -2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2;
     -2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2];
 
  a = 30;
  A11 = 30*eye(length(X));
  A12 = 0.0001*eye(length(X));
  A13 = [0.4 0 0 0;0 0.3 0 0;0 0 0.5 0];
 A3 = [A11;A12;A13];
 
% for j = 1:1:11
%     Ph(j) = exp(-norm(A11*X-C(:,j))^2/(2*b^2));
%       
% end
   
A112 = 0.02*rand(11,4);

% A112= [-4.9984   -4.9993   -4.9986   -4.9984;
%    -4.9998   -4.9991   -4.9988   -4.9984;
%    -4.9994   -4.9981   -5.0000   -4.9998;
%    -4.9980   -4.9988   -4.9980   -4.9997;
%    -4.9986   -4.9993   -4.9991   -4.9997;
%    -4.9981   -4.9983   -4.9983   -4.9997;
%    -4.9981   -4.9982   -4.9983   -4.9995;
%    -4.9986   -4.9986   -4.9992   -4.9981;
%    -4.9990   -4.9997   -4.9996   -4.9988;
%    -4.9986   -4.9998   -4.9983   -4.9984;
%    -4.9996   -4.9990   -4.9995   -4.9983];

Xa = [x11 x12 x21 x22]';
% Xa = [x11^2 x12^2 x21^2 x22^2 x11*x12 x11*x21 x11*x22 x12*x21 x12*x22 x21*x22 x11*x12*x21*x22]';
for j = 1:1:11
    V77 = rand(11,4);
    Ph(j) = 2/(1+exp(-V0(j,:)*X)) - 1;
%     Ph(j) = exp(-norm(A12*Xa-C(:,j))^2/(2*b^2));  
%     Ph(j) = (1-exp(-V0(j,:)*Xa/2))/(1 + exp(-V0(j,:)*Xa/2)) ;  
%       Ph(j) = (1 + cos(V0(j,:)*Xa))/2; 
end

% Ph=[x11^2 x12^2 x21^2 x22^2 x11*x12 x11*x21 x11*x22 x12*x21 x12*x22 x21*x22 x11*x12*x21*x22]';


% REFERENCE SIGNAL
Theta11d   = 0.4;
dTheta12d  = 0.2;
ddThetad = 0.1;

Psi11d   = 0.5;
dPsi12d  = 0.2;
ddPsi12d = 0.24;



Yd = [Theta11d dTheta12d Psi11d  dPsi12d]'; %reference vector
dYd = [dTheta12d ddThetad dPsi12d ddPsi12d]'; %derivative of reference vector

% f = [((k*r^2)/4*J1)*sin(x21);
%      ((k*r^2)/4*J2)*sin(x12)];


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
% 
% f = [(m1*g*r/J1-(k*r^2)/(J1))*sin(x11)+ ((k*r^2)/4*J1)*sin(x21) + ((k*r)/(2*J1))*(l-b);
%      (m2*g*r/J2-(k*r^2)/(J2))*sin(x21)+((k*r^2)/4*J2)*sin(x11) + ((k*r)/(2*J2))*(l-b)];



A = [0 1 0 0;
     0 0 0 0;        
     0 0 0 1;        
     0 0 0 0];

B = [0 0;
     1 0;
     0 0;
     0 1];
 
% fd  = [((k*r^2)/4*J1)*sin(y2d);
%        ((k*r^2)/4*J2)*sin(dy1d)];

% fd = [(m1*g*r/J1-(k*r^2)/(J1))*sin(y1d)+ ((k*r^2)/4*J1)*sin(y2d) + ((k*r)/(2*J1))*(l-b);
%       (m2*g*r/J2-(k*r^2)/(J2))*sin(y2d)+((k*r^2)/4*J2)*sin(y1d) + ((k*r)/(2*J2))*(l-b)];

Fz = B*f - (dYd - A*Yd);

% 
% L0 =60*eye(4);

L01 = 15; L02 = 15; L03=15; L04 = 15;
L0 = [L01 0 0 0;0 L02 0 0;0 0 L03 0;0 0 0 L04];

    
g0 =  [1/Iy     0;
       0     1/Iz]; 
  






Dhat = (L0)*(Z - Eta');

% if t==10
%     Dhat = (L0)*(Z - Eta')*exp(-0.5*sin(t));
% end


Hz = eye(4);

 No = [ 0.8 0.1 0 0;
        0.1 0.9 0 0;
         0  0 1 0;
         0  0 0 1];

%Calculation of sliding mode control
E1 = No*Hz*B*g0;
E2 = inv(E1'*E1)*E1';
alo = .10;
n1  = 5;
eta = 0.03;
R1 = S'*E1;
R11 = abs(R1(1));R12 = abs(R1(2));%Individual components of Right hand side of Eq. 20

Ks = [.107   .305    .208   .12;
      .105    .208    .302    .71];
  
% ud = inv(g0)*B'*( dYd - B*fd - A*Yd);
  
us = -W'*Ph- E2*No*Hz*Dhat - R*atan(E1'*S/eta)- Ks*S;


%% PI BASED DESIGN

Kp011 = 1.5; Kp012 = 1.2;
Kp021 = 1.2; Kp022 = 2; l0 = 4*[3.5 0;0 3.5]; l1 = .06*[0.3 0;0 0.3]; 

Kp0 = [Kp011 Kp012;Kp021 Kp022];
kep = 1;
Kp = l0*Theta - kep*Kp0;

alph = .5*[1 0 0 0;0 0 1 0];

KI = Kp*alph;

u0 = -(Kp0 + Kp)*Et;


% u0 = [0 0]';
%% Control Law
uc = u0 + us;
% ======================


% uc = ud + uz + us;

Gm1 = 1e-3*eye(11);
ko=2;
lm =0.7;
% Gm2 = 8e-5*rand*eye(10);
%Update Laws 
SEta = (A*Z + Fz + B*g0*uc + B*g0*W'*Ph + Dhat); %disturbance observer Eq.8
SR = alo*(R11+R12 - lm*R); %+ d0*norm(Z)+d1*R; %Disturbance gain update Eq. 20
SW = Gm1*(Ph*S'*E1/2 + n1*W);
 

SWG1 = SW(:,1); SWG2 = SW(:,2);
% 
% P0 = 0.05e-0;

GammaTheta = [0.02 0.01;0.01 0.02];

GTheta = GammaTheta*( l0*norm(Et)^2 * Kp0 - l1*Theta);

LL1 = 6e-1*eye(11);
Av = -10*eye(4);
rhv = 0.2;
nv = 1;

SV0 = LL1*((eye(11) - diag(V0*X))'*W * (B*g0)'* (inv(Av))'*Z *rhv * sign(X') + nv *norm(Z)*V0);%check here

for i=1:1:4
    sys(i) = SEta(i);
end

for i=5:1:5
    sys(i) = SR(i-4);
end

for i = 6:1:16  
    sys(i) = SWG1(i-5);
end

for i = 17:1:27  
    sys(i) = SWG2(i-16);
end

for i = 28:1:29
    sys(i) = GTheta(1,i-27);
end

for i = 30:1:31
    sys(i) = GTheta(2,i-29);
end

for i = 32:1:42
    sys(i) = SV0(i-31,1);
end

for i = 43:1:53
    sys(i) = SV0(i-42,2);
end

for i = 54:1:64
    sys(i) = SV0(i-53,3);
end

for i = 65:1:75
    sys(i) = SV0(i-64,4);
end



function sys=mdlOutputs(t,x,u)

global Ph

z1  = u(1);
z2  = u(2);
z3  = u(3);
z4  = u(4);
x11 = u(5);
x12 = u(6);
x21 = u(7);
x22 = u(8);
S1  = u(9);
S2  = u(10);
S3  = u(11);
S4  = u(12);
Et1 = u(13);
Et2 = u(14);
ml  = u(15);
hl  = u(16);
% Dz1 = u(15);
% Dz2 = u(16);
% Dz3 = u(17);
% Dz4 = u(18);
% 
% Dz = [Dz1 Dz2 Dz3 Dz4]';


Et = [Et1 Et2]';     
X = [x11 x12 x21 x22]';
S = [S1 S2 S3 S4]';



for i =1:1:4
    Eta(i)=x(i);
end

for i =1:1:1
    R(i) = x(4+i);
end

for i = 1:1:11  
    Wf(i) = x(5+i);
end

for i = 1:1:11
    Ws(i) = x(16+i);
end

W = [Wf' Ws'];

for i = 1:1:2
    Theta(1,i) = x(27+i);
end

for i = 1:1:2
    Theta(2,i) = x(29+i);
end
   
for i = 1:1:11
    V0(i,1) = x(31+i);
end

for i = 1:1:11
    V0(i,2) = x(42+i);
end

for i = 1:1:11
    V0(i,3) = x(53+i);
end

for i = 1:1:11
    V0(i,4) = x(64+i);
end



Z = [z1 z2 z3 z4]';


for j = 1:1:11
    V77 = rand(11,4);
    
    Ph(j) = 2/(1+exp(-V0(j,:)*X)) - 1;
%   Ph(j) = exp(-norm(A12*Xa-C(:,j))^2/(2*b^2));
%     Ph(j) = (1-exp(-V0(j,:)*Xa/2))/(1 + exp(-V0(j,:)*Xa/2));  
%       Ph(j) = (1 + cos(V0(j,:)*Xa))/2; 
end

% REFERENCE SIGNAL

% REFERENCE SIGNAL
Theta11d   = 0.4;
dTheta12d  = 0.2;
ddThetad = 0.1;

Psi11d   = 0.5;
dPsi12d  = 0.2;
ddPsi12d = 0.24;



Yd = [Theta11d dTheta12d Psi11d  dPsi12d]'; %reference vector
dYd = [dTheta12d ddThetad dPsi12d ddPsi12d]'; %derivative of reference vector

A = [0 1 0 0;
     0 0 0 0;        
     0 0 0 1;        
     0 0 0 0];

B = [0 0;
     1 0;
     0 0;
     0 1];
 

 
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







% f = [-(dml/M)*x12; PI1*x22];
% % 
      
fd = [PI2*dTheta12d; PI3*dPsi12d];

L01 = 15; L02 = 15; L03=15; L04 = 15;
L0 = [L01 0 0 0;0 L02 0 0;0 0 L03 0;0 0 0 L04];

g0 =  [1/Iy     0;
       0     1/Iz];  


Dhat = (L0)*(Z - Eta');

Hz = eye(4);
 
 No = [ 0.8 0.1 0 0;
        0.1 0.9 0 0;
         0  0 1 0;
         0  0 0 1];

%% Calculation of sliding mode control
E1 = No*Hz*B*g0;
E2 = inv(E1'*E1)*E1';
eta = 0.03;
Ks = [.107    .305    .208    .12;
      .105    .208    .302    .71];

us = -W'*Ph- E2*No*Hz*Dhat - R*atan(E1'*S/eta)- Ks*S;

ud = inv(g0)*B'*( dYd - B*fd - A*Yd); 


%% PI BASED DESIGN

Kp011 = 1.5; Kp012 = 1.2;
Kp021 = 1.2; Kp022 = 2; l0 = 4*[3.5 0;0 3.5]; 

Kp0 = [Kp011 Kp012;Kp021 Kp022];
kep = 1;
Kp = l0*Theta - kep*Kp0;

alph = .5*[1 0 0 0;0 0 1 0];

KI = Kp*alph;
u0 = -(Kp0 + Kp)*Et; %PI controller

%% Control Law
uc = u0 + us;
 %% Sloshing Forces
 %% The sloshing forces "Fs" are defined as follows: 
         %From equation 2.10, 2.11, and 2.12, we consider i=1, j=0 for the Fsx and 
         % i = 0 for Fsy. Then
omega_00 = 0.2; 
omega_10 = 0.4; 
Beta10 = 4; 
Beta00 = 2; 
lambda_10 = 0.2; 
lambda_00 = 0.3; 
rho_s = 1.002; 
Fsx = -rho_s * omega_10^2 *Beta10 * lambda_10 * sin(omega_10*t);
Fsy = -rho_s * omega_00^2 *Beta00 * lambda_00 * sin(omega_00*t);
Fs  = [Fsx Fsy 0]'; %Total sloshing forces


ax = 0.36;
K11 = 1;
Omega_sx = - rho_s * omega_10^2 * Beta10 * ax * sin(omega_10*t);      
Omega_sy = (rho_s * omega_00^2 * Beta00 * sin(omega_00*t))/K11;
Tau_s = [Omega_sx; Omega_sy; 0];
 
Delta = [Tau_s(2); Tau_s(3)];


sig1 = 0; sig2 = 0;

sigma = [sig1 0;
         0    sig2];

% d = [5*sin(0.65*pi*t+0.65)+5;
%      3.5*cos(0.63*pi*t+0.2)+3];



%% Esternal Disturbances

Va = 0.3;
Fh1 = 1.2; 
Fh2 = 2.3; 
Fh3 = 4.1; 
Fh4 = 1.6;%hub forces
b  = 0.3; %angle of sideslip of the quadrotor
hf = Fh1+Fh2+Fh3+Fh4; %hub forces generated by the 4 rotors
Fh = [hf*cos(b) hf*sin(b) 0]'; %total hub force vector

rho_a = 0.89; 
Sa = 0.4; 
Cxy = 0.2; 
Cz = 0.1;
Fd1 = -1/2 * rho_a * Sa*Cxy *cos(b)*Va^2;
Fd2 = -1/2 * rho_a * Sa*Cxy *sin(b)*Va^2;
Fd3 = 1/2 * rho_a * Sa*Cz*Va^2;
Fd = [Fd1 Fd2 Fd3]'; %Total air resistance vector

dr = Fh + Fd; %Total disturbance vector

%% The moment disturbances "dw = tau0" are defined as follows:
l = 0.225; %arm length
h = 0.2;
Tau_h1 = -h*(Fh1+Fh2+Fh3+Fh4)*sin(b);
Tau_h2 = h*(Fh1+Fh2+Fh3+Fh4)*cos(b);
Tau_h3 = l*(Fh2 - Fh4)*cos(b) + l*(Fh3 - Fh1)*sin(b);

Tau_h  = [Tau_h1 Tau_h2 Tau_h3]'; %Total torque due to hub forces

Jr     = 0.003; 
Omega_r = 0.03; 
p = 0.3; 
q = 0.4;
Tau_w1 = Jr*q*Omega_r;
Tau_w2 = -Jr*p*Omega_r;
Tau_w = [Tau_w1 Tau_w2 0]'; %Total torque due to the aerodynamic forces 
Tau_0 = Tau_h + Tau_w; %Total force torque

if t<=15
    d1 = Tau_0(2);
else 
    d1 = 1;
end

if t<=5
    d2 = Tau_0(3);
end

if t>5
    d2 = Tau_0(3);
end

if t>=15
    d2 = 1+0.3*exp(-2*t);
end
   
d = [d1 d2]';

%% ==========================================================

%% Integral squared errors

ISE1 = abs(sqrt((Z(1))^2+(Z(2))^2));%sum(Z(1) + Z(2))^2;
ISE2 = abs(sqrt((Z(3))^2+(Z(4))^2));%sum(Z(3) + Z(4))^2;

sys(1) = uc(1);
sys(2) = uc(2);
sys(3) = u0(1);
sys(4) = u0(2);
sys(5) = Kp(1,1);
sys(6) = Kp(1,2);
sys(7) = Kp(2,1);
sys(8) = Kp(2,2);
sys(9) = d(1);
sys(10)= Dhat(2);
sys(11)= KI(1,1);
sys(12)= KI(2,3);
% sys(11)= Dhat(4);
% sys(12)= d(2);
% sys(13)= KI(1,1);
% sys(14)= KI(1,2);
% sys(15)= KI(2,3);
% sys(16)= KI(2,4);
sys(13) = ISE1;
sys(14) = ISE2;
sys(15) = d(2);
sys(16)= Dhat(4);


