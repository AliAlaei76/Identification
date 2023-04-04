%% Deriving Cov(theta)

s=tf('s');
G=-6.2447/((s+0.2423)*(s^2+0.35*s+77.37));
t=0:0.1:20;
N=numel(t);
u=wgn(N,1,1);
[y1,t]=lsim(G,u,t);
y11=y1+0.005*rand(N,1); % Measured output
U=zeros(N,8);
y=[0;0;0;0;y11];
x=[0;0;0;u];
for i=1:N
    U(i,:)=[-y(i+3) -y(i+2) -y(i+1) -y(i) x(i+3) x(i+2) x(i+1) x(i)];
end
[P,R,Q1]=svd(U);
Q=Q1';
p=8;
ystar=P'*y11;
y1star=ystar(1:p);
Rbracket=R(1:p,1:p);
thetahatstar=Rbracket^(-1)*y1star;
thetahat=Q'*thetahatstar
yhat=U*thetahat
ehat=y11-yhat
sigmahat2=(1/(N-p))*ehat'*ehat
cov=sigmahat2*(U'*U)^(-1)