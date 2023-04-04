%% EMM method
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
p=21;
u1e=[0;0;0;0;0;u1];
y1e=[0;0;0;0;0;y1];
e=zeros(N+n,1);
v=zeros(N+n,1);
theta0=[-0.385977605433590
0.0710082831824715;
-0.555940425612344;
-0.0641808802380704;
-0.0277474946715707;
-0.0101064686383221;
0.445271787549417;
0.114300904599795;
-0.0580537563168732;
-0.189782687856995;
-0.327036905183232;
zeros(10,1)
];
theta=[theta0 zeros(p,N)];
%% Implementation of EMM method

ut=zeros(N,p);
P=zeros(p,p,N+1);
P(:,:,1)=0.1*eye(p);
for i=1:N
    ut(i,:)=[-y1e(i+4) -y1e(i+3) -y1e(i+2) -y1e(i+1) -y1e(i) u1e(i+5) u1e(i+4) u1e(i+3) u1e(i+2) u1e(i+1) u1e(i) -e(i+4) -e(i+3) -e(i+2) -e(i+1) -e(i) v(i+4) v(i+3) v(i+2) v(i+1) v(i)];
    P(:,:,i+1)=P(:,:,i)-(1/(1+ut(i,:)*P(:,:,i)*ut(i,:)'))*(P(:,:,i)*ut(i,:)'*ut(i,:)*P(:,:,i));
    theta(:,i+1)=theta(:,i)+P(:,:,i)*ut(i,:)'*(y1(i+5)-ut(i,:)*theta(:,i));
    u1t=ut(i,1:11);
    theta1=theta(1:11,i+1);
    e(i+5)=y1(i+5)-u1t*theta1;
    v(i+5)=y1(i+5)-ut(i,:)*theta(:,i+1);
    if norm(theta(:,i+1)-theta(:,i))<10^(-5)
        break
    end
end
numberofiteration=i+1
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

plot(m,theta(12,1:numberofiteration+1),m,theta(13,1:numberofiteration+1),m,theta(14,1:numberofiteration+1),m,theta(15,1:numberofiteration+1),m,theta(16,1:numberofiteration+1),'linewidth',1.5)
grid on
xlabel('Iteration')
ylabel('$\hat{\theta}$','Interpreter',"latex",'FontSize',13)
legend('d_1','d_2','d_3','d_4','d_5')

plot(m,theta(17,1:numberofiteration+1),m,theta(18,1:numberofiteration+1),m,theta(19,1:numberofiteration+1),m,theta(20,1:numberofiteration+1),m,theta(21,1:numberofiteration+1),'linewidth',1.5)
grid on
xlabel('Iteration')
ylabel('$\hat{\theta}$','Interpreter',"latex",'FontSize',13)
legend('f_1','f_2','f_3','f_4','f_5')