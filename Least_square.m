%% System identification based on least square

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
U
theta=(U.'*U)^(-1)*U.'*y11
a1=theta(1);
a2=theta(2);
a3=theta(3);
a4=theta(4);
b0=theta(5);
b1=theta(6);
b2=theta(7);
b3=theta(8);
%% Plot

z=tf('z',0.1);
G2=(b0+b1*z^(-1)+b2*z^(-2)+b3*z^(-3))/(1+a1*z^(-1)+a2*z^(-2)+a3*z^(-3)+a4*z^(-4))
[yhat,t]=lsim(G2,u,t);
plot(t,y11,t,yhat)
grid on
legend('Response of system','Response of identified model')
xlabel('Time(second)')
ylabel('y(t)')
%% Step plot

u2=ones(N,1);
[y2,t]=lsim(G,u2,t);
y22=y2+0.005*rand(N,1);
[yhat2,t]=lsim(G2,u2,t);
plot(t,y22,t,yhat2)
grid on
legend('step Response of system','step Response of identified model')
xlabel('Time(second)')
ylabel('S(t)')