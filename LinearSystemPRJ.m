clc;
clear all;
%% Variables
a = 1.65;
b = 6.2;
c = 2.69;
d = 4.2;
m = 0.12;
g = 9.8;
y1o = 0.02;
y2o = -0.02;
yc = 0.12;
y12o = yc + y2o - y1o;
u1o = a*((y1o + b)^4)*((c/(y12o + d)^4)+m*g);
u2o = a*((y2o + b)^4)*((-c/(y12o + d)^4)+m*g);
t=0:0.01:10;
%%
k1 = 4*u1o/(a*(y1o+b)^5);
k2 = 4*u2o/(a*(-y2o+b)^5);
k3 = 4*c/(y12o+d)^5;
k4 = 1/(a*(y1o+b)^4);
k5 = 1/(a*(-y2o+b)^4);
%% state space matrices
A = [0 1 0 0;
    -(k1+k3)/m 0 k3/m 0;
     0 0 0 1;
    k3/m 0 -(k2+k3)/m 0]
B = [0 0;
    k4/m 0;
     0 0;
     0 k5/m]
C = [1 0 0 0;
     0 0 1 0;]
D = zeros(size(C,1),size(B,2))
%%
sys = ss(A,B,C,D);
sys_tf = tf(sys)
%% Canonical Forms
% Observable Canonical Form
csys = canon(sys,'companion');
Ao=csys.A;
Bo=csys.B;
Co=csys.C;
Do=csys.D;
% Controllable Canonical Form
Ac=transpose(Ao);
Bc=transpose(Co);
Cc=transpose(Bo);
Dc=transpose(Do);
% Jordan Canonical Form
[V,J] = jordan(A);
%% Step Response & Impulse Response
t=0:0.01:10;
figure
step(sys_tf(1,1),t)
hold on;
figure
step(sys_tf(2,2),t)
hold on;
figure
impulse(sys_tf(1,1),t)
hold on
figure
impulse(sys_tf(2,2),t)
hold on
%% Checking Controllability and Obervability
sys_OL = sys;
sys_order = order(sys_OL)
ctrb_rank = rank(ctrb(A,B))
obsv_rank = rank(obsv(A,C))
%% Bode Diagram
figure
bode(sys_OL(1,1))
hold on;
figure
bode(sys_OL(2,2))
hold on;
%% Root locus
figure
rlocus(sys(1,1))
hold on;
figure
rlocus(sys(2,2))
hold on;
%% Check open loop eigenvalues
E=eig(A)
%Solve for K using pole placement
P=[-2.4+3.2i -2.4-3.2i -2.4+3.2i -2.4-3.2i];    %Put 4 poles in desired place
Kc = place(A,B,P);
%Create closed loop system
sys_cl = ss(A-B*Kc,B,C,D);
%Check closed loop eigenvalues
EE=eig(A-B*Kc)
%Check step response
t = 0:0.01:5;
figure
step(sys_cl,t)
grid
title('Step Response with State-Feedback Controller')
%% Estimator
Pe = 4*[-2.4+3.2i -2.4-3.2i -2+3i -2-3i];
obsv(A,C)
g1 = 53248  - 39.4817;
g2 = 8089.6 - 0;
g3 = 771.2  - 12.5699;
g4 = 35.2  - 0;
G = [g1 g1;g2 g2;g3 g3;g4 g4];
s=tf('s');
Ges = Kc*(s*eye(4)-A+B*Kc+G*C)*G
%% Reduced order estimator
A1e = A(1,2:end);
Ae1 = A(2:end,1);
Aee = A(2:end,2:end);
b1 = B(1,:);
Be = B(2:end,:);
