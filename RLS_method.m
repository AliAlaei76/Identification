%% RLS method
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
U=zeros(20,11);
for i=1:20
    U(i,:)=[-y1e(i+4) -y1e(i+3) -y1e(i+2) -y1e(i+1) -y1e(i) u1e(i+5) u1e(i+4) u1e(i+3) u1e(i+2) u1e(i+1) u1e(i)];
end
theta0=(U'*U)^(-1)*U'*y1(1:20);
P0=100*eye(p);
%% Implemention of RLS method

shift=20;
theta=[theta0 zeros(p,N)];
P=zeros(p,p,N);
P(:,:,1)=P0;
ut=zeros(N,p);
for i=21:N
    ut(i-shift,:)=[-y1(i-1) -y1(i-2) -y1(i-3) -y1(i-4) -y1(i-5) u1(i) u1(i-1) u1(i-2) u1(i-3) u1(i-4) u1(i-5)];
    P(:,:,i-shift+1)=P(:,:,i-shift)-(1/(1+ut(i-shift,:)*P(:,:,i-shift)*ut(i-shift,:)'))*(P(:,:,i-shift)*ut(i-shift,:)'*ut(i-shift,:)*P(:,:,i-shift));
    theta(:,i-shift+1)=theta(:,i-shift)+P(:,:,i-shift+1)*ut(i-shift,:)'*(y1(i)-ut(i-shift,:)*theta(:,i-shift));
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

yhat=zeros(i,1);
for j=21:i
    yhat(j)=ut(j-shift,:)*theta(:,j-shift+1);
end
e=y1(21:N)-yhat(21:N);
s=e'*e
sigma=s/(N-p)

plot(t(21:N),y1(21:N),t(21:N),yhat(21:N))
grid on
xlabel('Time(second)')
ylabel('y(t)')
legend('y(t)','yhat(t)')

plot(t(21:N),e)
grid on
xlabel('Time(second)','FontSize',15);
ylabel('$\hat{e}$','Interpreter',"latex",'FontSize',15)