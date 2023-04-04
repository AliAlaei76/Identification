%% Deriving system order

clear;
clc;
close all;
s=tf('s');
G=-6.2447/((s+0.2423)*(s^2+0.35*s+77.37));
t=0:0.1:20;
N=numel(t);
u1=wgn(N,1,1);
[y,t]=lsim(G,u1,t);
y1=y+0.005*rand(N,1); % Measured output
%% First order

N=numel(y1);
p11=3;
u11=[0;u1];
y11=[0;y1];
U=zeros(N,3);
for i=1:N
    U(i,:)=[-y11(i) u11(i+1) u11(i)];
end
theta11=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y11hat=U*theta11;
e11=y1-y11hat;
s11=e11'*e11
sigma11=s11/(N-p11)
%% Second order

p12=5;
u12=[0;0;u1];
y12=[0;0;y1];
U=zeros(N,5);
for i=1:N
    U(i,:)=[-y12(i+1) -y12(i) u12(i+2) u12(i+1) u12(i)];
end
theta12=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y12hat=U*theta12;
e12=y1-y12hat;
s12=e12'*e12
sigma12=s12/(N-p12)
%% Third order

p13=7;
u13=[0;0;0;u1];
y13=[0;0;0;y1];
U=zeros(N,7);
for i=1:N
    U(i,:)=[-y13(i+2) -y13(i+1) -y13(i) u13(i+3) u13(i+2) u13(i+1) u13(i)];
end
theta13=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y13hat=U*theta13;
e13=y1-y13hat;
s13=e13'*e13
sigma13=s13/(N-p13)
%% Fourth order

p14=9;
u14=[0;0;0;0;u1];
y14=[0;0;0;0;y1];
U=zeros(N,9);
for i=1:N
    U(i,:)=[-y14(i+3) -y14(i+2) -y14(i+1) -y14(i) u14(i+4) u14(i+3) u14(i+2) u14(i+1) u14(i)];
end
theta14=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y14hat=U*theta14;
e14=y1-y14hat;
s14=e14'*e14
sigma14=s14/(N-p14)
%% Fifth order

p15=11;
u15=[0;0;0;0;0;u1];
y15=[0;0;0;0;0;y1];
U=zeros(N,11);
for i=1:N
    U(i,:)=[-y15(i+4) -y15(i+3) -y15(i+2) -y15(i+1) -y15(i) u15(i+5) u15(i+4) u15(i+3) u15(i+2) u15(i+1) u15(i)];
end
theta15=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y15hat=U*theta15;
e15=y1-y15hat;
s15=e15'*e15
sigma15=s15/(N-p15)
%% Sixth order

p16=13;
u16=[0;0;0;0;0;0;u1];
y16=[0;0;0;0;0;0;y1];
U=zeros(N,13);
for i=1:N
    U(i,:)=[-y16(i+5) -y16(i+4) -y16(i+3) -y16(i+2) -y16(i+1) -y16(i) u16(i+6) u16(i+5) u16(i+4) u16(i+3) u16(i+2) u16(i+1) u16(i)];
end
theta16=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y16hat=U*theta16;
e16=y1-y16hat;
s16=e16'*e16
sigma16=s16/(N-p16)
%% Seventh order

p17=15;
u17=[0;0;0;0;0;0;0;u1];
y17=[0;0;0;0;0;0;0;y1];
U=zeros(N,15);
for i=1:N
    U(i,:)=[-y17(i+6) -y17(i+5) -y17(i+4) -y17(i+3) -y17(i+2) -y17(i+1) -y17(i) u17(i+7) u17(i+6) u17(i+5) u17(i+4) u17(i+3) u17(i+2) u17(i+1) u17(i)];
end
theta17=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y17hat=U*theta17;
e17=y1-y17hat;
s17=e17'*e17
sigma17=s17/(N-p17)
%% Eighth order

p18=17;
u18=[0;0;0;0;0;0;0;0;u1];
y18=[0;0;0;0;0;0;0;0;y1];
U=zeros(N,17);
for i=1:N
    U(i,:)=[-y18(i+7) -y18(i+6) -y18(i+5) -y18(i+4) -y18(i+3) -y18(i+2) -y18(i+1) -y18(i) u18(i+8) u18(i+7) u18(i+6) u18(i+5) u18(i+4) u18(i+3) u18(i+2) u18(i+1) u18(i)];
end
theta18=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y18hat=U*theta18;
e18=y1-y18hat;
s18=e18'*e18
sigma18=s18/(N-p18)
%% Ninth order

p19=19;
u19=[0;0;0;0;0;0;0;0;0;u1];
y19=[0;0;0;0;0;0;0;0;0;y1];
U=zeros(N,19);
for i=1:N
    U(i,:)=[-y19(i+8) -y19(i+7) -y19(i+6) -y19(i+5) -y19(i+4) -y19(i+3) -y19(i+2) -y19(i+1) -y19(i) u19(i+9) u19(i+8) u19(i+7) u19(i+6) u19(i+5) u19(i+4) u19(i+3) u19(i+2) u19(i+1) u19(i)];
end
theta19=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y19hat=U*theta19;
e19=y1-y19hat;
s19=e19'*e19
sigma19=s19/(N-p19)
%% Tenth order

p10=21;
u10=[0;0;0;0;0;0;0;0;0;0;u1];
y10=[0;0;0;0;0;0;0;0;0;0;y1];
U=zeros(N,21);
for i=1:N
    U(i,:)=[-y10(i+9) -y10(i+8) -y10(i+7) -y10(i+6) -y10(i+5) -y10(i+4) -y10(i+3) -y10(i+2) -y10(i+1) -y10(i) u10(i+10) u10(i+9) u10(i+8) u10(i+7) u10(i+6) u10(i+5) u10(i+4) u10(i+3) u10(i+2) u10(i+1) u10(i)];
end
theta10=(U'*U)^(-1)*U'*y1;
deg=rank(U'*U)
y10hat=U*theta10;
e10=y1-y10hat;
s10=e10'*e10
sigma10=s10/(N-p10)
%% Plot

n=1:1:10;
s=[s11 s12 s13 s14 s15 s16 s17 s18 s19 s10]
bar(n,s,0.4)
grid on
xlabel('n','FontSize',15);
ylabel('$\hat{s}$','Interpreter',"latex",'FontSize',15)
sigma=[sigma11 sigma12 sigma13 sigma14 sigma15 sigma16 sigma17 sigma18 sigma19 sigma10]
bar(n,sigma,0.4)
grid on
xlabel('n','FontSize',15);
ylabel('$\hat{\sigma}^2$','Interpreter',"latex",'FontSize',15)
plot(t,e11,t,e15)
grid on
xlabel('Time(second)','FontSize',15);
ylabel('$\hat{e}$','Interpreter',"latex",'FontSize',15)
legend('error of first order estimation','error of Fifth order estimation')