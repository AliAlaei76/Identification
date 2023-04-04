%% Continuous state space representation of system
A=[0 1 0;-77.3707 -0.2703 -6.2447;-0.9879 0 -0.322];
B=[0 0 1].';
C=[1 0 0];
D=0;
sys=ss(A,B,C,D); x0=[0.5468 0 9.6926].';
%% Continuous transfer function
s=tf('s');
G=C*(s*eye(3)-A)^(-1)*B+D
p=pole(G)
zpk(sys)
%% Continuous impulse response
syms s;
G1=-6.245/(s^3+0.5923*s^2+77.46*s+18.74)
g=ilaplace(G1,s)
g=vpa(g,2)
impulse(G)
grid on
%% Continuous frequency response
bode(sys)
grid on
%% Discrete state space
dss=c2d(sys,1)
step(dss)
grid on
%% Discrete transfer function
Gd=c2d(G,1)
%% discrete impulse response
syms z
Gd1=(-0.06697*z^2-0.1053*z-0.04708)/(z^3+0.5714*z^2-0.3598*z-0.5531);
h=iztrans(Gd1)
h=vpa(h,2)
impulse(Gd)
%% discrete frequency response
bode(dss)
grid on