%% System identification based on correlation

s=tf('s');
G=-6.2447/((s+0.2423)*(s^2+0.35*s+77.37));
t=0:0.1:20;
N=numel(t);
u=wgn(N,1,1);
[y1,t]=lsim(G,u,t);
y11=y1+0.005*rand(N,1); % Measured output
%% Deriving Rxx , Rxy

y=[y11;y11];
Rxy=zeros(N,1);
sum=0;
for k=1:N
    for i=1:N
        sum=sum+u(i)*y(i+(k-1));
    end
    Rxy(k)=(1/N)*sum;
    sum=0;
end

Rxx=zeros(N,1);
sum=0;
x=[u;u];
for k=1:N
    for i=1:N
        sum=sum+x(i)*x(i+(k-1));
    end
    Rxx(k)=(1/N)*sum;
    sum=0;
end
%% Deriving H(z)

z=tf('z',0.1);
rxyz=0;
for k=1:N
    rxyz=Rxy(k)*z^(-(k-1))+rxyz;
end
rxxz=0;
for k=1:N
    rxxz=Rxx(k)*z^(-(k-1))+rxxz;
end
H=rxyz/rxxz;
%% Plot

[yhat,t]=lsim(H,u,t);
plot(t,y11,t,yhat)
grid on
legend('Response of system','Response of identified model')
xlabel('Time(second)')
ylabel('y(t)')