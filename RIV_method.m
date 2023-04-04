%% RIV method
%% Estimate of Theta0 and P0

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
n=5;
p=11;
u1e=[0;0;0;0;0;u1];
y1e=[0;0;0;0;0;y1];
U=zeros(N,11);
for i=1:N
    U(i,:)=[-y1e(i+4) -y1e(i+3) -y1e(i+2) -y1e(i+1) -y1e(i) u1e(i+5) u1e(i+4) u1e(i+3) u1e(i+2) u1e(i+1) u1e(i)];
end
theta0=(U(1:20,:)'*U(1:20,:))^(-1)*U(1:20,:)'*y1(1:20);
P0=100*eye(p);
shift=20;
yhat=U(1:20,:)*theta0;
yhat=[yhat;zeros(N-shift,1)];
P0=100*eye(p);
%% Implemention of RIV method

shift=20;
theta=[theta0 zeros(p,N)];
P=zeros(p,p,N);
P(:,:,1)=P0;
Z=zeros(N,p);
ut=zeros(N,p);
for i=21:N
    ut(i-shift,:)=[-y1(i-1) -y1(i-2) -y1(i-3) -y1(i-4) -y1(i-5) u1(i) u1(i-1) u1(i-2) u1(i-3) u1(i-4) u1(i-5)];
    Z(i-shift,:)=[-yhat(i-1) -yhat(i-2) -yhat(i-3) -yhat(i-4) -yhat(i-5) u1(i) u1(i-1) u1(i-2) u1(i-3) u1(i-4) u1(i-5)];
    P(:,:,i-shift+1)=P(:,:,i-shift)-(1/(1+ut(i-shift,:)*P(:,:,i-shift)*Z(i-shift,:)'))*(P(:,:,i-shift)*Z(i-shift,:)'*ut(i-shift,:)*P(:,:,i-shift));
    theta(:,i-shift+1)=theta(:,i-shift)+P(:,:,i-shift+1)*Z(i-shift,:)'*(y1(i)-ut(i-shift,:)*theta(:,i-shift));
    yhat(i)=ut(i-shift,:)*theta(:,i-shift+1);
    if norm(theta(:,i-shift+1)-theta(:,i-shift))<10^(-6)
        break
    end
end
numberofiteration=i-shift
%% Plots and results

m=0:1:numberofiteration;
plot(m,theta(1,1:numberofiteration+1),m,theta(2,1:numberofiteration+1),m,theta(3,1:numberofiteration+1),m,theta(4,1:numberofiteration+1),m,theta(5,1:numberofiteration+1),'linewidth',1.5)
grid on
xlabel('Iteration')
ylabel('$\hat{\theta}$','Interpreter',"latex",'FontSize',13)
legend('a_1','a_2','a_3','a_4','a_5')

plot(m,theta(6,1:numberofiteration+1),m,theta(7,1:numberofiteration+1),m,theta(8,1:numberofiteration+1),m,theta(9,1:numberofiteration+1),m,theta(10,1:numberofiteration+1),m,theta(11,1:numberofiteration+1),'linewidth',1.5)
grid on
xlabel('Iteration')
ylabel('$\hat{\theta}$','Interpreter',"latex",'FontSize',13)
legend('b_0','b_1','b_2','b_3','b_4','b_5')

ehat=y1-yhat;
shat=ehat'*ehat
sigma=shat/(N-p)

plot(t,y1,t,yhat)
grid on
xlabel('Time(second)')
ylabel('y(t)')
legend('y(t)','yhat(t)')

plot(t,ehat)
grid on
xlabel('Time(second)','FontSize',15);
ylabel('$\hat{e}$','Interpreter',"latex",'FontSize',15)