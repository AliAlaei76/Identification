%% System identification via least square method based on impulse response

s=tf('s');
G=-6.2447/((s+0.2423)*(s^2+0.35*s+77.37));
t=0:0.1:20;
N=numel(t);
u=wgn(N,1,1);
[y1,t]=lsim(G,u,t);
y11=y1+0.005*rand(N,1); % Measured output
M=100;
U=zeros(N,M);
x=[zeros((M-1),1);u];
for i=1:N
    for j=1:M
        U(i,j)=x(i+M-j);
    end
end
theta=(U.'*U)^(-1)*U.'*y11
%% Plot

u2=[1;zeros(M-1,1)];
[y2,t(1:M)]=lsim(G,u2,t(1:M));
plot(t(1:M),y2,t(1:M),theta)
grid on
legend('impulse response of system','impulse response of model')
xlabel('time(second)')
ylabel('h(t)')