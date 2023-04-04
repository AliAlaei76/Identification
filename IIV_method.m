%% IIV method

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
theta0=(U'*U)^(-1)*U'*y1
iter=30;
theta=[theta0 zeros(p,iter)];
yhat=zeros(N,iter);
yhate=zeros(N+n,iter);
Z=zeros(N,p,iter);
for k=1:iter
    yhat(:,k)=U*theta(:,k);
    yhate(n+1:N+n,k)=yhat(:,k);
    for i=1:N
        Z(i,:,k)=[-yhate(i+4,k) -yhate(i+3,k) -yhate(i+2,k) -yhate(i+1,k) -yhate(i,k) u1e(i+5) u1e(i+4) u1e(i+3) u1e(i+2) u1e(i+1) u1e(i)];
    end
    theta(:,k+1)=(Z(:,:,k)'*U)^(-1)*Z(:,:,k)'*y1;
end
%% Plots and results

m=0:1:iter;
plot(m,theta(1,:),m,theta(2,:),m,theta(3,:),m,theta(4,:),m,theta(5,:),'linewidth',1.5)
grid on
xlabel('Iteration')
ylabel('$\hat{\theta}$','Interpreter',"latex",'FontSize',12)
legend('a_1','a_2','a_3','a_4','a_5')
plot(m,theta(6,:),m,theta(7,:),m,theta(8,:),m,theta(9,:),m,theta(10,:),m,theta(11,:),'linewidth',1.5)
grid on
xlabel('Iteration')
ylabel('$\hat{\theta}$','Interpreter',"latex",'FontSize',13)
legend('b_0','b_1','b_2','b_3','b_4','b_5')

yhatlast=U*theta(:,iter+1);
ehat=y1-yhatlast;
shat=ehat'*ehat
sigma=shat/(N-p)
plot(t,y1,t,yhatlast)
grid on
xlabel('Time(second)')
ylabel('y(t)')
legend('y(t)','yhat(t)')
plot(t,ehat)
grid on
xlabel('Time(second)','FontSize',15);
ylabel('$\hat{e}$','Interpreter',"latex",'FontSize',15)