function [sys,x0,str,ts] = Control_UpdateLaws(t,x,u,flag) 
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
sizes.NumContStates = 247;
sizes.NumDiscStates = 0;
sizes.NumOutputs = 25; 
sizes.NumInputs = 12; 
sizes.DirFeedthrough = 2; 
sizes.NumSampleTimes = 0; 
sys = simsizes(sizes); 


y1d0   = sin(0)+sin(0.5*(0));
dy1d0  = cos(0)+(0.5)*cos(0.5*(0));
y2d0   = cos(0)+cos(0);
dy2d0  = -sin(0)-sin(0);

X0  =[0.6 -0.6 2 -2]';



% Kp011 = 0.4; Kp012 = 0.5;
% Kp021 = 0.2; Kp022 = 0.7;



Yd0 = [y1d0 dy1d0 y2d0 dy2d0]';
Eta = X0 -Yd0;

% Eta = [4 -3 -2 2]';%Initializing Eta of disturbance observer in Eq. 8

R = .10;%initializing the disturbance compensation gain of Eq. 20

% Wf = 5e-3*ones(11,1);
% Ws = 5e-3*ones(11,1);
% Wf=-1+(1+1)*rand(11,1);
% Wf  = 5e-1*rand(11,1);%first weight for RBF (11 by 1) of Eq. 23
% Ws  = 3*rand(11,1);%second weight for RBF (11 by 1) of Eq.23

Wf=   [0.9018
   -0.1121
   -0.8800
    0.7335
    0.2624
   -0.2899
    0.9940
   -0.5517
    0.3049
    0.2100
   -0.2255]';
Ws  = Wf;%second weight for RBF (11 by 1) of Eq.23

% Wf  = 0.0005*rand(11,1);%first weight for RBF (11 by 1) of Eq. 23
% Ws  = 0.003*rand(11,1);%second weight for RBF (11 by 1) of Eq.23
Ph  = [0 0 0 0 0 0 0 0 0 0 0]';%Initializing RBF (11 by 1) of Eq. 23

A1 = zeros(10,10);
% A1 = [0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 0 0 0;
%       0 0 0 0 0 0 0 0 0 0;];%initializing Lambda1 of Eq. 64

A2 = zeros(1,10);
           
% A2 = [0 0 0 0 0 0 0 0 0 0];%Initializing Lambda2 of Eq. 65

% W1 = 4-2*rand(10,1); %weight for ADP (10 by 1) of Eq. 71
% %  W1 = .05e-1*ones(10,1);
%  W1=-1+(1+1)*rand(10,1);
%  W1=-1+(-1+5)*rand(10,1);

W1 =[4.2163
    4.5070
    0.9799
   -0.3738
    2.9734
    0.4432
    0.7925
    2.8745
    4.9460
    0.6390]';

% W1 =[4.2163
%     4.5070
%     0.9799
%    -0.3738
%     2.9734
%     0.4432
%     0.7925
%     2.8745
%     4.9460
%     0.6390]';%works great, random between -1 and 5;

% W1 =[-0.6060
%     1.4265
%     1.3929
%     0.2134
%     1.0007
%    -0.6453
%     0.0544
%     1.2893
%     0.9805
%     1.3987]';%2nd Best
% 
% W1 =[0.9843
%     1.7036
%     0.8862
%    -0.0298
%     1.9482
%     0.1265
%     0.4779
%     2.0480
%     2.2713
%     1.6524]';%It works really well

% W1=[-0.0183 -0.0221 -0.1854 -0.18854 0.5559 ... 
%      0.6559 1.2316 0.0212 0.0212 0.1851]';
Gamma = eye(10);%initializing Gamma of Eq. 71

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

for iij = 1:1:10
    x0(27+iij) = A1(1,iij);
end

for iiij = 1:1:10
    x0(37+iiij) = A1(2,iiij);
end

for ijj = 1:1:10
    x0(47+ijj) = A1(3,ijj);
end

for ijjj = 1:1:10
    x0(57+ijjj) = A1(4,ijjj);
end

for jj = 1:1:10
    x0(67+jj) = A1(5,jj);
end

for jjj = 1:1:10
    x0(77+jjj) = A1(6,jjj);
end

for ji = 1:1:10
    x0(87+ji) = A1(7,ji);
end

for jii = 1:1:10
    x0(97+jii) = A1(8,jii);
end

for jji = 1:1:10
    x0(107+jji) = A1(9,jji);
end

for jjji = 1:1:10
    x0(117+jjji) = A1(10,jjji);
end

for k = 1:1:10
    x0(127+k) = A2(k);
end

for kk = 1:1:10
    x0(137+kk) = W1(kk);
end

for ki = 1:1:10
    x0(147+ki) = Gamma(1,ki);
end

for ik = 1:1:10
    x0(157+ik) = Gamma(2,ik);
end

for kii = 1:1:10
    x0(167+kii) = Gamma(3,kii);
end

for kik = 1:1:10
    x0(177+kik) = Gamma(4,kik);
end

for iki = 1:1:10
    x0(187+iki) = Gamma(5,iki);
end

for iik = 1:1:10
    x0(197+iik) = Gamma(6,iik);
end

for kiik = 1:1:10
    x0(207+kiik) = Gamma(7,kiik);
end

for id = 1:1:10
    x0(217+id) = Gamma(8,id);
end

for di = 1:1:10
    x0(227+di) = Gamma(9,di);
end

for ddi = 1:1:10
    x0(237+ddi) = Gamma(10,ddi);
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

for ie = 1:1:10
   A1(1,ie) = x(27+ie);
end

for iee = 1:1:10
   A1(2,iee) = x(37+iee);
end


for ei = 1:1:10
    A1(3,ei) = x(47+ei);
end

for iv = 1:1:10
    A1(4,iv) = x(57+iv);
end


for vi = 1:1:10
    A1(5,vi) = x(67+vi);
end

for iq = 1:1:10
    A1(6,iq) = x(77+iq);
end

for qi = 1:1:10
    A1(7,qi) = x(87+qi);
end

for ti = 1:1:10
    A1(8,ti) = x(97+ti);
end

for hi = 1:1:10
    A1(9,hi) = x(107+hi);
end

for ih = 1:1:10
    A1(10,ih) = x(117+ih);
end

L1 = [A1(1,:);A1(2,:);A1(3,:);A1(4,:);A1(5,:);
       A1(6,:);A1(7,:);A1(8,:);A1(9,:);A1(10,:)];

for iy = 1:1:10
    A2(iy) = x(127+iy);
end
L2 = A2;

for iyy = 1:1:10
    W1(iyy) = x(137+iyy);
end
W1 =W1';

for yyi = 1:1:10
    Gamma(1,yyi) = x(147+yyi);
end

for tti = 1:1:10
    Gamma(2,tti) = x(157+tti);
end

for lli = 1:1:10
    Gamma(3,lli) =x(167+lli);
end

for i = 1:1:10
    Gamma(4,i) = x(177+i);
end

for i = 1:1:10
    Gamma(5,i) = x(187+i);
end

for i = 1:1:10
    Gamma(6,i) = x(197+i);
end

for i = 1:1:10
    Gamma(7,i) = x(207+i);
end

for i = 1:1:10
    Gamma(8,i) = x(217+i);
end

for i = 1:1:10
    Gamma(9,i) = x(227+i);
end

for i = 1:1:10
    Gamma(10,i) = x(237+i);
end

Gm = 0.1e-3*[Gamma(1,:);Gamma(2,:);Gamma(3,:);Gamma(4,:);Gamma(5,:);
      Gamma(6,:);Gamma(7,:);Gamma(8,:);Gamma(9,:);Gamma(10,:)];

   
Z = [z1 z2 z3 z4]';

b = 2;

% C = -2+(2+2)*rand(4,11);
C = [-2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2;
     -2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2;
     -2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2;
     -2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2];
 
  a = 30;
  A11 = 0.30*eye(length(X));
  A12 = 0.3*eye(length(X));
  A13 = [0.4 0 0 0;0 0.3 0 0;0 0 0.5 0];
 A3 = [A11;A12;A13];
 
%  Ph2 = [x11^2 x12^2 x21^2 x22^2 x11*x12 x11*x21 x11*x22 x12*x21 x12*x22 x21*x22 x11*x12*x21*x22];
 Ph2 = [z1^2 z2^2 z3^2 z4^2 z1*z2 z1*z3 z1*z4 z2*z3 z2*z4 z3*z4 z1*z2*z3*z4];
for j = 1:1:11
%     Ph(j) = exp(-norm(A11*X-C(:,j))^2/(2*b^2));
%       Ph(j) = Ph2(j);
      Ph(j) = 2/(1+exp(-A3(j,:)*X)) - 1;
      
end

% Ph1 = [z1^2 z2^2 z3^2 z4^2 z1*z2 z1*z3 z1*z4 z2*z3 z2*z4 z3*z4]';%Activation function
% 
%    
% dPh1 = [2*z1   0  0    0;
%         0   2*z2  0    0;
%         0    0   2*z3  0;
%         0    0    0    2*z4;
%         z2   z1   0    0;
%         z3   0    z1   0;
%         z4   0    0    z1;
%         0    z3   z2   0;
%         0    z4   0    z2;
%         0    0    z4   z3]; %Derivative of activation func wrt z1,z2,z3,z4


Ph1 = [z1^2 z1*z2 z1*z3 z1*z4 z2^2 ...
       z2*z3 z2*z4 z3^2 z3*z4 z4^2]'; %Activation function

dPh1 = [2*z1   0  0  0;
        z2   z1  0   0;
        z3   0   z1  0;
        z4   0   0   z1;
        0   2*z2 0   0;
        0    z3  z2  0;
        0    z4  0   z2;
        0    0  2*z3 0;
        0    0   z4  z3;
        0    0   0   2*z4]; %Derivative of activation func wrt z1,z2,z3,z4
    
    

% REFERENCE SIGNAL

y1d   = sin(t)+sin(0.5*t);
dy1d  = cos(t)+(0.5)*cos(0.5*t);
ddy1d = -sin(t)-(0.5)*(0.5)*sin(0.5*t);

y2d   = cos(t)+cos(t);
dy2d  = -sin(t)-sin(t);
ddy2d = -cos(t)-cos(t);

Yd = [y1d dy1d y2d dy2d]'; %reference vector
d0 = norm(Yd);%bound on desired signal 
d1 = 0.002; %constant term to attenuate the disturbance

dYd = [dy1d ddy1d dy2d ddy2d]'; %derivative of reference vector

k = 10; J1 = 1; J2 = 1; m1 =2; m2 = 2; r = 0.1;l = 0.5; b = 0.4; g = 9.8;

f = [(m1*g*r/J1-(k*r^2)/(J1))*sin(x11)+ ((k*r^2)/4*J1)*sin(x21) + ((k*r)/(2*J1))*(l-b);
     (m2*g*r/J2-(k*r^2)/(J2))*sin(x21)+((k*r^2)/4*J2)*sin(x11) + ((k*r)/(2*J2))*(l-b)];

A = [0 1 0 0;
     0 0 0 0;        
     0 0 0 1;        
     0 0 0 0];

B = [0 0;
     1 0;
     0 0;
     0 1];
 
fd = [(m1*g*r/J1-(k*r^2)/(J1))*sin(y1d)+ ((k*r^2)/4*J1)*sin(y2d) + ((k*r)/(2*J1))*(l-b);
      (m2*g*r/J2-(k*r^2)/(J2))*sin(y2d)+((k*r^2)/4*J2)*sin(y1d) + ((k*r)/(2*J2))*(l-b)];
   
Fz = B*f - (dYd - A*Yd);
  
Fd = B*(f - fd)+(B*B' - eye(4))*( dYd - A*Yd);

% ADP controller 
% Ro = [ .01 0;
%        0 .01];
Ro = [0.2 0;0 0.2]; 


% L0 =  [10 0 0 0;
%        0 10 0 0;
%        0 0 12 0;
%        0 0 0 10];

% L0 =60*eye(4);
L01 = 100; L02 = 100; L03=100; L04 = 100;
L0 = [L01 0 0 0;0 L02 0 0;0 0 L03 0;0 0 0 L04];

% Gm1 = diag(10,10);
% Gm1 = 0.0080*rand*eye(11);

% Gm1 =[10 0 0 0 0 0 0 0 0 0 0;
%       0 11 0 0 0 0 0 0 0 0 0;
%       0 0 12 0 0 0 0 0 0 0 0;
%       0 0 0 13 0 0 0 0 0 0 0;
%       0 0 0 0 15 0 0 0 0 0 0;
%       0 0 0 0 0 16 0 0 0 0 0;
%       0 0 0 0 0 0 17 0 0 0 0;
%       0 0 0 0 0 0 0 18 0 0 0;
%       0 0 0 0 0 0 0 0 19 0 0;
%       0 0 0 0 0 0 0 0 0 20 0;
%       0 0 0 0 0 0 0 0 0 0 21];
    
g0 = [1/J1 0;
      0 1/J2];
  
%Desired control
ud = inv(g0)*B'*(dYd - B*fd - A*Yd);

%n3 = sin(t)^2*cos(t)+sin(t)^5 + sin(1.12*t)

% if t<=4
%     uz = -(1/2)*inv(Ro)*g0'*B'*dPh1'*W1 + [sin(t)*cos(0.02*t) sin(.4*t)*cos(0.02*t)]';
%     
% else
     uz = -(1/2)*inv(Ro)*g0'*B'*dPh1'*W1;
% end

Jhat = W1'*Ph1;




%Weight update law
Y = dPh1*(A*Z +Fd +(1/2)*B*g0*uz);
Yo = dPh1*(A*Z +Fd +B*g0*uz);
% Q = [ .20 .8 0 0;
%       .8 .30 0 0;
%       0  0 1 0;
%       0  0 0 1];

Q = .1*eye(4);

Xo = Z'*Q*Z;

Dhat = L0*(Z - Eta');
al2 = 3e-2;
Beta2 = 10;

Hz = eye(4);
 
 No = [ 0.8 0.1 0 0;
        0.1 0.9 0 0;
         0  0 1 0;
         0  0 0 1];

%Calculation of sliding mode control
E1 = No*Hz*B*g0;
E2 = inv(E1'*E1)*E1';
alo = 20;
n1  = .1;
eta = 0.03;
R1 = S'*E1;
R11 = abs(R1(1));R12 = abs(R1(2));%Individual components of Right hand side of Eq. 20

Ks = [0.7250    0.5308    0.8706    0.2991;
      0.5872    0.8154    0.2823    0.1279];
  
us = -W'*Ph- E2*No*Hz*Dhat - R*atan(E1'*S/eta)- Ks*S;

%Control Law

% uc = ud+uz+us;
uc = ud + uz + us;

Gm1 = 8e-1*eye(11);
ko=2;
% Gm2 = 8e-5*rand*eye(10);
%Update Laws 
SEta = (A*Z + Fz + B*g0*uc + B*g0*W'*Ph + Dhat); %disturbance observer Eq. 8
SR = alo*(R11+R12); %+ d0*norm(Z)+d1*R; %Disturbance gain update Eq. 20
SW = Gm1*(Ph*S'*E1 + n1*norm(Z/ko)*W);

% To = Y*Y';
% if min(eig(To))<0
%     To = To+rand()*eye(length(Y));
%     SL1  = -al2*L1 + To;
%     SL2  = -al2*L2 + Xo*Y';
% else
    SL1  = -al2*L1 + Y*Y';
    SL2  = -al2*L2 + Xo*Y';
% end


% SL1  = -al2*L1 + Y*Y';
% SL2  = -al2*L2 + Xo*Y';
Qo = 2*Z;
Go = B*g0*inv(Ro)*g0'*B';

% if Qo'*(A*Z+Fd+B*g0*uz)<=0
%     Sw1  = -Gm*( L2' + L1'*W1);
% %     Sw1  =-al2*Yo/((1+Yo'*Yo)^2)*(Xo+W1'*dPh1*(A*Z + Fd)-(1/4)*W1'*dPh1*B*g0*inv(Ro)*g0'*B'*dPh1'*W1);
% else
    

    Pi =1; yo=0.5;
%     Sw1  =-al2*Yo/((1+Yo'*Yo)^2)*(Xo+W1'*dPh1*(A*Z + Fd)-(1/4)*W1'*dPh1*B*g0*inv(Ro)*g0'*B'*dPh1'*W1)+yo/2 * Pi*dPh1*Go*Qo;
%     Sw1  =-al2*Yo/((1+Yo'*Yo)^2)*(Xo+W1'*dPh1*(A*Z + Fd)-(1/4)*W1'*dPh1*B*g0*inv(Ro)*g0'*B'*dPh1'*W1);
%     Sw1  = -Gm*( L2' + L1'*W1)+yo/2 * (Pi*dPh1*Go*Qo);
% end
%     
    
% Sw1  = -Gm2*( L2' - L1'*W1)+yo/2 * Pi*Go*Qo;%
Sw1  = -al2*Yo/((1+Yo'*Yo)^2)*(Xo+W1'*dPh1*(A*Z + Fd)-(1/4)*W1'*dPh1*B*g0*inv(Ro)*g0'*B'*dPh1'*W1);

Gmo = [4.0e-11,0,0,0,0,0,0,0,0,0;
       0,4.0e-11,0,0,0,0,0,0,0,0;
       0,0,4.0e-11,0,0,0,0,0,0,0;
       0,0,0,4.0e-11,0,0,0,0,0,0;
       0,0,0,0,4.0e-11,0,0,0,0,0;
       0,0,0,0,0,4.0e-11,0,0,0,0;
       0,0,0,0,0,0,4.0e-11,0,0,0;
       0,0,0,0,0,0,0,4.0e-11,0,0;
       0,0,0,0,0,0,0,0,4.0e-11,0;
       0,0,0,0,0,0,0,0,0,4.0e-11];
   
   ro = 5;
   ro1 = 0.04;

% if min(eig(L1))<=ro
%     SGm =  - Gm*L1*Gm;
% else
    SGm  = Beta2*Gm - ro1*Gm*(Y*Y'/ro)*Gm;
% end
 

SWG1 = SW(:,1); SWG2 = SW(:,2);

SL1G1 = SL1(1,:);SL1G2 = SL1(2,:);SL1G3 = SL1(3,:);SL1G4 = SL1(4,:);
SL1G5 = SL1(5,:);SL1G6 = SL1(6,:);SL1G7 = SL1(7,:);SL1G8 = SL1(8,:);
SL1G9 = SL1(9,:);SL1G10 = SL1(10,:);

GGm1 = SGm(1,:);GGm2 = SGm(2,:); GGm3 = SGm(3,:);GGm4 = SGm(4,:);GGm5 = SGm(5,:);
GGm6 = SGm(6,:);GGm7 = SGm(7,:); GGm8 = SGm(8,:);GGm9 = SGm(9,:);GGm10 = SGm(10,:);

% assert(all(imag(SL1)==0),'* is imaginary or nan')
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

for i = 28:1:37
   sys(i) = SL1G1(i-27);
end

for i = 38:1:47
   sys(i) = SL1G2(i-37);
end

for i = 48:1:57
   sys(i) = SL1G3(i-47);
end

for i = 58:1:67
   sys(i) = SL1G4(i-57);
end

for i = 68:1:77
   sys(i) = SL1G5(i-67);
end

for i = 78:1:87
   sys(i) = SL1G6(i-77);
end

for i = 88:1:97
   sys(i) = SL1G7(i-87);
end

for i = 98:1:107
   sys(i) = SL1G8(i-97);
end

for i = 108:1:117
   sys(i) = SL1G9(i-107);
end

for i = 118:1:127
   sys(i) = SL1G10(i-117);
end

for i = 128:1:137
    sys(i) = SL2(i-127);
end

for i = 138:1:147
    sys(i) = Sw1(i-137);
end

for i = 148:1:157
    sys(i) = GGm1(i-147);
end

for i = 158:1:167
    sys(i) = GGm2(i-157);
end

for i = 168:1:177
    sys(i) = GGm3(i-167);
end

for i = 178:1:187
    sys(i) = GGm4(i-177);
end

for i = 188:1:197
    sys(i) = GGm5(i-187);
end

for i = 198:1:207
    sys(i) = GGm6(i-197);
end

for i = 208:1:217
    sys(i) = GGm7(i-207);
end

for i = 218:1:227
    sys(i) = GGm8(i-217);
end

for i = 228:1:237
    sys(i) = GGm9(i-227);
end

for i = 238:1:247
    sys(i) = GGm10(i-237);
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


Delta = [2*cos( x11 + x21 )^2 + 30*x21;
         0.4*cos( x11 + x21 )^2 + ( (2*x11^2 + 4*x21^3)/(1+(x11+x21)^2 ) ) + exp(-x21^2)];

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

for i = 1:1:10
   A1(1,i) = x(27+i);
end

for i = 1:1:10
   A1(2,i) = x(37+i);
end


for i = 1:1:10
    A1(3,i) = x(47+i);
end

for i = 1:1:10
    A1(4,i) = x(57+i);
end


for i = 1:1:10
    A1(5,i) = x(67+i);
end

for i = 1:1:10
    A1(6,i) = x(77+i);
end

for i = 1:1:10
    A1(7,i) = x(87+i);
end

for i = 1:1:10
    A1(8,i) = x(97+i);
end

for i = 1:1:10
    A1(9,i) = x(107+i);
end

for i = 1:1:10
    A1(10,i) = x(117+i);
end

L1 = [A1(1,:);A1(2,:);A1(3,:);A1(4,:);A1(5,:);
       A1(6,:);A1(7,:);A1(8,:);A1(9,:);A1(10,:)];

for i = 1:1:10
    A2(i) = x(127+i);
end
L2 = A2;

for i = 1:1:10
    W1(i) = x(137+i);
end
W1 =W1';

% for i = 1:1:10
%     Gamma(1,i) = x(147+i);
% end
% 
% for i = 1:1:10
%     Gamma(2,i) = x(157+i);
% end
% 
% for i = 1:1:10
%     Gamma(3,i) =x(167+i);
% end
% 
% for i = 1:1:10
%     Gamma(4,i) = x(177+i);
% end
% 
% for i = 1:1:10
%     Gamma(5,i) = x(187+i);
% end
% 
% for i = 1:1:10
%     Gamma(6,i) = x(197+i);
% end
% 
% for i = 1:1:10
%     Gamma(7,i) = x(207+i);
% end
% 
% for i = 1:1:10
%     Gamma(8,i) = x(217+i);
% end
% 
% for i = 1:1:10
%     Gamma(9,i) = x(227+i);
% end
% 
% for i = 1:1:10
%     Gamma(10,i) = x(237+i);
% end

% Gm = [Gamma(1,:);Gamma(2,:);Gamma(3,:);Gamma(4,:);Gamma(5,:);
%       Gamma(6,:);Gamma(7,:);Gamma(8,:);Gamma(9,:);Gamma(10,:)];

   
Z = [z1 z2 z3 z4]';

b = 2;
% C = -2+(2+2)*rand(4,11);

C = [-2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2;
     -2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2;
     -2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2;
     -2 -1.6 -1.2 -0.8 -0.4 0 0.4 0.8 1.2 1.6 2];
 
  a = 30;
  A11 = 0.30*eye(length(X));
  A12 = 0.3*eye(length(X));
  A13 = [0.4 0 0 0;0 0.3 0 0;0 0 0.5 0];
 A3 = [A11;A12;A13];
 
%  Ph2 = [x11^2 x12^2 x21^2 x22^2 x11*x12 x11*x21 x11*x22 x12*x21 x12*x22 x21*x22 x11*x12*x21*x22];
 Ph2 = [z1^2 z2^2 z3^2 z4^2 z1*z2 z1*z3 z1*z4 z2*z3 z2*z4 z3*z4 z1*z2*z3*z4];
for j = 1:1:11
%     Ph(j) = exp(-norm(A11*X-C(:,j))^2/(2*b^2));
%       Ph(j) = Ph2(j);
      Ph(j) = 2/(1+exp(-A3(j,:)*X)) - 1;
      
end


Deltahat = W'*Ph; 


Delta1Tilder = Delta(1)-Deltahat(1);
Delta2Tilder = Delta(2)-Deltahat(2);
% Ph1 = [z1^2 z2^2 z3^2 z4^2 z1*z2 z1*z3 z1*z4 z2*z3 z2*z4 z3*z4]';%Activation function
% 
%    
% dPh1 = [2*z1   0  0    0;
%         0   2*z2  0    0;
%         0    0   2*z3  0;
%         0    0    0    2*z4;
%         z2   z1   0    0;
%         z3   0    z1   0;
%         z4   0    0    z1;
%         0    z3   z2   0;
%         0    z4   0    z2;
%         0    0    z4   z3]; %Derivative of activation func wrt z1,z2,z3,z4


Ph1 = [z1^2 z1*z2 z1*z3 z1*z4 z2^2 ...
       z2*z3 z2*z4 z3^2 z3*z4 z4^2]'; %Activation function

dPh1 = [2*z1   0  0  0;
        z2   z1  0   0;
        z3   0   z1  0;
        z4   0   0   z1;
        0   2*z2 0   0;
        0    z3  z2  0;
        0    z4  0   z2;
        0    0  2*z3 0;
        0    0   z4  z3;
        0    0   0   2*z4]; %Derivative of activation func wrt z1,z2,z3,z4
    

% REFERENCE SIGNAL

y1d   = sin(t)+sin(0.5*t);
dy1d  = cos(t)+(0.5)*cos(0.5*t);
ddy1d = -sin(t)-(0.5)*(0.5)*sin(0.5*t);

y2d   = cos(t)+cos(t);
dy2d  = -sin(t)-sin(t);
ddy2d = -cos(t)-cos(t);

Yd = [y1d dy1d y2d dy2d]'; %reference vector
dYd = [dy1d ddy1d dy2d ddy2d]'; %derivative of reference vector

k = 10; J1 = 1; J2 = 1; m1 =2; m2 = 2; r = 0.1;l = 0.5; b = 0.4; g = 9.8;

A = [0 1 0 0;
     0 0 0 0;        
     0 0 0 1;        
     0 0 0 0];

B = [0 0;
     1 0;
     0 0;
     0 1];
 
fd = [(m1*g*r/J1-(k*r^2)/(J1))*sin(y1d)+ ((k*r^2)/4*J1)*sin(y2d) + ((k*r)/(2*J1))*(l-b);
      (m2*g*r/J2-(k*r^2)/(J2))*sin(y2d)+((k*r^2)/4*J2)*sin(y1d) + ((k*r)/(2*J2))*(l-b)];

% ADP controller 
% Ro = [ .01 0;
%        0 .01];

Ro = [0.2 0;0 0.2]; 
% Ro = [0.2 0;0 0.2];
% L0 =  [10 0 0 0;
%        0 10 0 0;
%        0 0 12 0;
%        0 0 0 10];

% L0 =60*eye(4);
L01 = 100; L02 = 100; L03=100; L04 = 100;
L0 = [L01 0 0 0;0 L02 0 0;0 0 L03 0;0 0 0 L04];

g0 = [1/J1 0;
      0 1/J2];    
%Desired control

ud = inv(g0)*B'*( dYd - B*fd - A*Yd);    
    
% if t<=4
%     uz = -(1/2)*inv(Ro)*g0'*B'*dPh1'*W1 + [sin(t)*cos(0.2*t) sin(.4*t)*cos(0.2*t)]';
%     
% else
     uz = -(1/2)*inv(Ro)*g0'*B'*dPh1'*W1;
% end

Jhat = W1'*Ph1;

%Weight update law

Dhat = L0*(Z - Eta');

Hz = eye(4);
 
 No = [ 0.8 0.1 0 0;
        0.1 0.9 0 0;
         0  0 1 0;
         0  0 0 1];

%Calculation of sliding mode control
E1 = No*Hz*B*g0;
E2 = inv(E1'*E1)*E1';
eta = 0.03;
Ks = [0.7250    0.5308    0.8706    0.2991;
      0.5872    0.8154    0.2823    0.1279];

us = -W'*Ph- E2*No*Hz*Dhat - R*atan(E1'*S/eta)- Ks*S;



%% PI based design



%Control Law
% uc = ud+uz+us;
uc = ud + uz + us;

% Del = -uc + ud + uz - R*sign(E1'*S) - E2*No*Hz*Dhat - Ks*S;
DDhat = pinv((E2*No*Hz)'*(E2*No*Hz))*(E2*No*Hz)'*( -uc +ud + uz - W'*Ph - R*atan(E1'*S/eta));

Del = -uc + ud + uz - R*sign(E1'*S) - E2*No*Hz*DDhat - Ks*S;


k = 10; J1 = 1; J2 = 1;
m1 =2; m2 = 2; r = 0.1;l = 0.5; b = 0.4; g = 9.8;

% d = [(m1*g*r/J1-(k*r^2)/(4*J1))*sin(x11)+((k*r)/(2*J1))*(l-b);
%      (m2*g*r/J2-(k*r^2)/(4*J2))*sin(x21)+((k*r)/(2*J2))*(l-b)];

% d = [5*sin(0.65*pi*t+0.65)+5;
%      3.5*cos(0.63*pi*t+0.2)+3];

if t<=15
    d1 = 5*sin(2*t + 1);
else 
    d1 = 1;
end

if t<=5
    d2 = 0;
end

if t>5
    d2 = 0.5*sin( 2*t );
end

if t>=15
    d2 = 1+0.3*exp(-2*t);
end
   
d = [d1 d2]';


Delta = [130*cos(x21 + x22)^2 + 3*x12;
         120*cos(x21 + x22)^2 + ((x21^4 + x12^3)/((5+x11+x12)^2 )) + exp(-x12^2)];

% GH = W'*Ph;     
% DH = Delta - Del;
%  B1 = normalize(W1);

%% Error performance
ISE1a = abs(sqrt((Z(1))^2+(Z(2))^2));%sum(Z(1) + Z(2))^2;
ISE2a = abs(sqrt((Z(3))^2+(Z(4))^2));%sum(Z(3) + Z(4))^2;



ed1 = d(1)-Dhat(2);
 
Wscale = W1 - min(W1);
EE = max(W1) - min(W1);
B1 = Wscale/EE;
lambda_min = min(eig(L1)); 
sys(1)=uc(1);
sys(2)=uc(2);
sys(3)=ud(1);
sys(4)=ud(2);
sys(5)=uz(1);
sys(6)=uz(2);
sys(7)=lambda_min;
sys(8)=d(1);
sys(9)=W1(1);
sys(10)=W1(2);
sys(11)=W1(3);
sys(12)=W1(4);
sys(13)=W1(5);
sys(14)=W1(6);
sys(15)=W1(7);
sys(16)=W1(8);
sys(17)=W1(9);
sys(18)=W1(10);
sys(19)=Dhat(4);
sys(20)=d(2);
sys(21)=Del(1);
sys(22)=Delta(1);
sys(23)=Jhat;
sys(24) = ISE1a;
sys(25) = ISE2a;


