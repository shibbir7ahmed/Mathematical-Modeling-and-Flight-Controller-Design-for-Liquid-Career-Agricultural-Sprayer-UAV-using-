function [sys,x0,str,ts] = Control_Laws(t,x,u,flag) 
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
sizes.NumOutputs = 19; 
sizes.NumInputs = 12; 
sizes.DirFeedthrough = 2; 
sizes.NumSampleTimes = 0; 
sys = simsizes(sizes);

y1d0   = sin(0)+sin(0.5*(0));
dy1d0  = cos(0)+(0.5)*cos(0.5*(0));
y2d0   = cos(0)+cos(0);
dy2d0  = -sin(0)-sin(0);

X0  =[0.6 -0.6 2 -2]';


Yd0 = [y1d0 dy1d0 y2d0 dy2d0]';
Eta = X0 -Yd0;

R = 0;%initializing the disturbance compensation gain of Eq. 20
% [-0.0183 -0.0221 -0.1854 -0.18854 -0.1954 0.5559 0.6559 1.2316 0.0212 0.0212 0.1851];

Wf  = 5e-2*rand(11,1);%first weight for RBF (11 by 1) of Eq. 23
Ws  = 0.3*rand(11,1);%second weight for RBF (11 by 1) of Eq.23
% Wf=[-0.0183 -0.0221 -0.1854 -0.18854 -0.1954 0.5559 ... 
%      0.6559 1.2316 0.0212 0.0212 0.1851]';
% Ws=[-0.0183 -0.0221 -0.1854 -0.18854 -0.1954 0.5559 ... 
%      0.6559 1.2316 0.0212 0.0212 0.1851]';
% Wf  = -5+(5+5)*rand(11,1);%first weight for RBF (11 by 1) of Eq. 23
% Ws  = -5+(5+5)*rand(11,1);%second weight for RBF (11 by 1) of Eq.23
Ph  = [0 0 0 0 0 0 0 0 0 0 0]';%Initializing RBF (11 by 1) of Eq. 23

% W1  = 0.04*rand(10,1); %weight for ADP (10 by 1) of Eq. 71
% W1 = (10,1);
% W1=[-0.0183 -0.0221 -0.1854 -0.18854 0.5559 ... 
%      0.6559 1.2316 0.0212 0.0212 0.1851]';
%  W1 = rand(10,1);

%  W1=-1+(1+1)*rand(10,1);
W1 =[0.9843
    1.7036
    0.8862
   -0.0298
    1.9482
    0.1265
    0.4779
    2.0480
    2.2713
    1.6524]';

A1 = zeros(10,10);%Initializing Lambda 1

A2 = zeros(1,10);%Initializing Lamb

Gamma = 5*eye(10);%initializing Gamma of Eq. 71
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

for kk = 1:1:10
    x0(27+kk) = W1(kk);
end

for iij = 1:1:10
    x0(37+iij) = A1(1,iij);
end

for iiij = 1:1:10
    x0(47+iiij) = A1(2,iiij);
end

for ijj = 1:1:10
    x0(57+ijj) = A1(3,ijj);
end

for ijjj = 1:1:10
    x0(67+ijjj) = A1(4,ijjj);
end

for jj = 1:1:10
    x0(77+jj) = A1(5,jj);
end

for jjj = 1:1:10
    x0(87+jjj) = A1(6,jjj);
end

for ji = 1:1:10
    x0(97+ji) = A1(7,ji);
end

for jii = 1:1:10
    x0(107+jii) = A1(8,jii);
end

for jji = 1:1:10
    x0(117+jji) = A1(9,jji);
end

for jjji = 1:1:10
    x0(127+jjji) = A1(10,jjji);
end

for k = 1:1:10
    x0(137+k) = A2(k);
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

x11  = u(1);
x12  = u(2);
x21  = u(3);
x22  = u(4);
z1   = u(5);
z2   = u(6);
z3   = u(7);
z4   = u(8);
S1   = u(9);
S2   = u(10);
S3   = u(11);
S4   = u(12);

Z = [z1 z2 z3 z4]';
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

for iyy = 1:1:10
    W1(iyy) = x(27+iyy);
end
W1 =W1';

for ie = 1:1:10
   A1(1,ie) = x(37+ie);
end

for iee = 1:1:10
   A1(2,iee) = x(47+iee);
end

for ei = 1:1:10
    A1(3,ei) = x(57+ei);
end

for iv = 1:1:10
    A1(4,iv) = x(67+iv);
end

for vi = 1:1:10
    A1(5,vi) = x(77+vi);
end

for iq = 1:1:10
    A1(6,iq) = x(87+iq);
end

for qi = 1:1:10
    A1(7,qi) = x(97+qi);
end

for ti = 1:1:10
    A1(8,ti) = x(107+ti);
end

for hi = 1:1:10
    A1(9,hi) = x(117+hi);
end

for ih = 1:1:10
    A1(10,ih) = x(127+ih);
end

L1 = [A1(1,:);A1(2,:);A1(3,:);A1(4,:);A1(5,:);
       A1(6,:);A1(7,:);A1(8,:);A1(9,:);A1(10,:)];
   
for iy = 1:1:10
    A2(iy) = x(137+iy);
end

L2 = A2;

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

Gm = 0.1e-2*[Gamma(1,:);Gamma(2,:);Gamma(3,:);Gamma(4,:);Gamma(5,:);
      Gamma(6,:);Gamma(7,:);Gamma(8,:);Gamma(9,:);Gamma(10,:)];

b = 2;

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
    Ph(j) = exp(-norm(A11*X-C(:,j))^2/(2*b^2));
%       Ph(j) = Ph2(j);
%       Ph(j) = 2/(1+exp(-A3(j,:)*X)) - 1;
      
end
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
    
% Ro = [ 1 0;
%        0 1];
 
Ro = 1*eye(2);

y1d   = sin(t)+sin(0.5*t);
dy1d  = cos(t)+(0.5)*cos(0.5*t);
ddy1d = -sin(t)-(0.5)*(0.5)*sin(0.5*t);

y2d   = cos(t)+cos(t);
dy2d  = -sin(t)-sin(t);
ddy2d = -cos(t)-cos(t);

Yd = [y1d dy1d y2d dy2d]'; %reference vector
dYd = [dy1d ddy1d dy2d ddy2d]'; %derivative of reference vector


k = 10; J1 = 1; J2 = 1; m1 =2; m2 = 2; r = 0.1;l = 0.5; b = 0.4; g = 9.8;
g0 = [1/J1 0;
      0 1/J2];
  
  
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

f = [(m1*g*r/J1-(k*r^2)/(J1))*sin(x11)+ ((k*r^2)/4*J1)*sin(x21) + ((k*r)/(2*J1))*(l-b);
     (m2*g*r/J2-(k*r^2)/(J2))*sin(x21)+((k*r^2)/4*J2)*sin(x11) + ((k*r)/(2*J2))*(l-b)];
   
 
  Fz = B*f - (dYd - A*Yd);

Fd = B*(f - fd)+(B*B' - eye(4))*( dYd - A*Yd);
  
% L0 =  1e-2*[2 0 0 0;
%             0 1 0 0;
%             0  0 1 0;
%             0  0 0 1];

L0 = 110*eye(4);
 
uz = -(1/2)*inv(Ro)*g0'*B'*dPh1'*W1;   
%Weight update law
Y = dPh1*(A*Z +Fd +(1/2)*B*g0*uz);
Yo = dPh1*(A*Z +Fd +B*g0*uz);
% Q = [ 2 0.8 0 0;
%       0.8 3 0 0;
%       0  0 1 0;
%       0  0 0 1];

Q = eye(4);
  
Xo = Z'*Q*Z;

Dhat = L0*(Z - Eta');
al2 = 3e-1;
Beta2 = 10;

% al2 = 0.3;
% Beta2 = 50;

 Hz = eye(4);
 
 No = [ 0.8 0.1 0 0;
        0.1 0.9 0 0;
         0  0 1 0;
         0  0 0 1];
 
 
%Calculation of sliding mode control
E1 = No*Hz*B*g0;
E2 = inv(E1'*E1)*E1';
alo = 5;
n1  = 1;
eta = 0.03;
R1 = S'*E1;
R11 = abs(R1(1));R12 = abs(R1(2));%Individual components of Right hand side of Eq. 20
Ks = [0.7250    0.5308    0.8706    0.2991;
      0.5872    0.8154    0.2823    0.1279];
  
us = -W'*Ph - R*atan(E1'*S/eta) - E2*No*Hz*Dhat - Ks*S;

%Desired control
ud = inv(g0)*B'*(dYd - B*fd - A*Yd);


Jhat = W1'*Ph1;

uc = us + ud + uz;


Gm1 = 8e-3*eye(11);
ko=0.2;
% Gm2 = 8e-5*rand*eye(10);
%Update Laws 
SEta = A*Z + Fz + B*g0*uc + B*g0*W'*Ph + Dhat; %disturbance observer Eq. 8
SR = alo*(R11+R12); %+ d0*norm(Z)+d1*R; %Disturbance gain update Eq. 20
SW = Gm1*(Ph*S'*E1 + n1*norm(Z)*W);

SWG1 = SW(:,1); SWG2 = SW(:,2);

Qo = 2*Z;
Go = B*g0*inv(Ro)*g0'*B';

% if Qo'*(A*Z+Fd+B*g0*uz)<=0
%     Sw1  = -Gm*( L2' + L1'*W1);
% %     Sw1  =-al2*Yo/((1+Yo'*Yo)^2)*(Xo+W1'*dPh1*(A*Z + Fd)-(1/4)*W1'*dPh1*B*g0*inv(Ro)*g0'*B'*dPh1'*W1);
% else
    

    Pi =1; yo=0.5;
%     Sw1  =-al2*Yo/((1+Yo'*Yo)^2)*(Xo+W1'*dPh1*(A*Z + Fd)-(1/4)*W1'*dPh1*B*g0*inv(Ro)*g0'*B'*dPh1'*W1)+yo/2 * Pi*dPh1*Go*Qo;
    Sw1  = -Gm*( L2' + L1'*W1)+yo/2 * (Pi*dPh1*Go*Qo);
% end
%     


SL1  = -al2*L1 + Y*Y';
SL2  = -al2*L2 + Xo*Y';


   ro = 5;
   ro1 = 0.4;

% if min(eig(L1))<=ro
%     SGm =  - Gm*L1*Gm;
% else
    SGm  = Beta2*Gm - ro1*Gm*(Y*Y'/ro)*Gm;
% end

SL1G1 = SL1(1,:);SL1G2 = SL1(2,:);SL1G3 = SL1(3,:);SL1G4 = SL1(4,:);
SL1G5 = SL1(5,:);SL1G6 = SL1(6,:);SL1G7 = SL1(7,:);SL1G8 = SL1(8,:);
SL1G9 = SL1(9,:);SL1G10 = SL1(10,:);

GGm1 = SGm(1,:);GGm2 = SGm(2,:); GGm3 = SGm(3,:);GGm4 = SGm(4,:);GGm5 = SGm(5,:);
GGm6 = SGm(6,:);GGm7 = SGm(7,:); GGm8 = SGm(8,:);GGm9 = SGm(9,:);GGm10 = SGm(10,:);





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
    sys(i) = Sw1(i-27);
end

for i = 38:1:47
   sys(i) = SL1G1(i-37);
end

for i = 48:1:57
   sys(i) = SL1G2(i-47);
end

for i = 58:1:67
    sys(i) = SL1G3(i-57);
end

for i = 68:1:77
   sys(i) = SL1G4(i-67);
end

for i = 78:1:87
   sys(i) = SL1G5(i-77);
end

for i = 88:1:97
   sys(i) = SL1G6(i-87);
end

for i = 98:1:107
   sys(i) = SL1G7(i-97);
end

for i = 108:1:117
   sys(i) = SL1G8(i-107);
end

for i = 118:1:127
   sys(i) = SL1G9(i-117);
end

for i = 128:1:137
   sys(i) = SL1G10(i-127);
end

for i = 138:1:147
    sys(i) = SL2(i-137);
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

x11  = u(1);
x12  = u(2);
x21  = u(3);
x22  = u(4);
z1   = u(5);
z2   = u(6);
z3   = u(7);
z4   = u(8);
S1   = u(9);
S2   = u(10);
S3   = u(11);
S4   = u(12);

X = [x11 x12 x21 x22]';
Z = [z1 z2 z3 z4]';
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

for iyy = 1:1:10
    W1(iyy) = x(27+iyy);
end
W1 =W1';


for ie = 1:1:10
   A1(1,ie) = x(37+ie);
end

for iee = 1:1:10
   A1(2,iee) = x(47+iee);
end

for ei = 1:1:10
    A1(3,ei) = x(57+ei);
end

for iv = 1:1:10
    A1(4,iv) = x(67+iv);
end

for vi = 1:1:10
    A1(5,vi) = x(77+vi);
end

for iq = 1:1:10
    A1(6,iq) = x(87+iq);
end

for qi = 1:1:10
    A1(7,qi) = x(97+qi);
end

for ti = 1:1:10
    A1(8,ti) = x(107+ti);
end

for hi = 1:1:10
    A1(9,hi) = x(117+hi);
end

for ih = 1:1:10
    A1(10,ih) = x(127+ih);
end

L1 = [A1(1,:);A1(2,:);A1(3,:);A1(4,:);A1(5,:);
       A1(6,:);A1(7,:);A1(8,:);A1(9,:);A1(10,:)];
   
for iy = 1:1:10
    A2(iy) = x(137+iy);
end

L2 = A2;

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

Gm = 8e-8*[Gamma(1,:);Gamma(2,:);Gamma(3,:);Gamma(4,:);Gamma(5,:);
      Gamma(6,:);Gamma(7,:);Gamma(8,:);Gamma(9,:);Gamma(10,:)];


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
    Ph(j) = exp(-norm(A11*X-C(:,j))^2/(2*b^2));
%       Ph(j) = Ph2(j);
%       Ph(j) = 2/(1+exp(-A3(j,:)*X)) - 1;
      
end


% Ph1 = [z1^2 z1*z2 z1*z3 z1*z4 z2^2 ...
%        z2*z3 z2*z4 z3^2 z3*z4 z4^2]'; %Activation function

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
    



Ro = 1*eye(2);

y1d   = sin(t)+sin(0.5*t);
dy1d  = cos(t)+(0.5)*cos(0.5*t);
ddy1d = -sin(t)-(0.5)*(0.5)*sin(0.5*t);

y2d   = cos(t)+cos(t);
dy2d  = -sin(t)-sin(t);
ddy2d = -cos(t)-cos(t);

Yd = [y1d dy1d y2d dy2d]'; %reference vector
dYd = [dy1d ddy1d dy2d ddy2d]'; %derivative of reference vector


k = 10; J1 = 1; J2 = 1; m1 =2; m2 = 2; r = 0.1;l = 0.5; b = 0.4; g = 9.8;

g0 = [1/J1 0;
      0 1/J2];
  
  
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

f = [(m1*g*r/J1-(k*r^2)/(J1))*sin(x11)+ ((k*r^2)/4*J1)*sin(x21) + ((k*r)/(2*J1))*(l-b);
     (m2*g*r/J2-(k*r^2)/(J2))*sin(x21)+((k*r^2)/4*J2)*sin(x11) + ((k*r)/(2*J2))*(l-b)];
   
 
Fz = B*f - (dYd - A*Yd);


% L0 =  1e-2*[2 0 0 0;
%             0 1 0 0;
%             0 0 1 0;
%             0 0 0 1];
L0 = 110*eye(4);
    
%Weight update law
% Y = dPh1*(A*Z +Fd +(1/2)*B*g0*uz);
% Q = [ 2 .8 0 0;
%       .8 3 0 0;
%       0  0 1 0;
%       0  0 0 1];
% Xo = Z'*Q*Z;

Dhat = L0*(Z - Eta');

% al2 = 0.3;
% Beta2 = 50;

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


%Desired control
ud = inv(g0)*B'*(dYd - B*fd - A*Yd);


uz = -(1/2)*inv(Ro)*g0'*B'*dPh1'*W1;

Jhat = W1'*Ph1;

uc = us + ud + uz;

DDhat = pinv((E2*No*Hz)'*(E2*No*Hz))*(E2*No*Hz)'*( -uc +ud + uz - W'*Ph - R*atan(E1'*S/eta));


m1 =2; m2 = 2;l = 0.5;g = 9.8;

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
 
 Delta = [2*cos( x11 + x21 )^2 + 60*x21;
         60*cos( x11 + x21 )^2 +  (2*x11^2 + 4*x21^3)/(1+(x11+x21)^2 )  + exp(-x21^2)];

 Deltahat = W'*Ph; 

sys(1)=uc(1);
sys(2)=uc(2);
sys(3)=R;
sys(4)=uz(1);
sys(5)=uz(2);
sys(6)=us(1);
sys(7)=us(2);
sys(8)=W1(1);
sys(9)=W1(2);
sys(10)=W1(3);
sys(11)=W1(4);
sys(12)=W1(5);
sys(13)=W1(6);
sys(14)=W1(7);
sys(15)=W1(8);
sys(16)=W1(9);
sys(17)=W1(10);
sys(18)=Dhat(4);
sys(19)=d(2);



