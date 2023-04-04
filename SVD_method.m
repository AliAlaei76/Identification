%% Deriving unknown parameters based on SVD method

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
y1star=ystar(1:p)
Rbracket=R(1:p,1:p)
thetahatstar=Rbracket^(-1)*y1star
thetahat=Q'*thetahatstar
a1=thetahat(1);
a2=thetahat(2);
a3=thetahat(3);
a4=thetahat(4);
b0=thetahat(5);
b1=thetahat(6);
b2=thetahat(7);
b3=thetahat(8);
%% Step plot

z=tf('z',0.1);
G2=(b0+b1*z^(-1)+b2*z^(-2)+b3*z^(-3))/(1+a1*z^(-1)+a2*z^(-2)+a3*z^(-3)+a4*z^(-4))
u2=ones(N,1);
[y2,t]=lsim(G,u2,t);
y22=y2+0.005*rand(N,1);
[yhat2,t]=lsim(G2,u2,t);
plot(t,y22,t,yhat2)
grid on
legend('step Response of system','step Response of identified model')
xlabel('Time(second)')
ylabel('S(t)')