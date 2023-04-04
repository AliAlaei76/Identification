%% System identification based on prony method

s=tf('s');
G=-6.2447/((s+0.2423)*(s^2+0.35*s+77.37));
t=0:1:40;
N=numel(t);
u=[1;zeros(N-1,1)];
[y,t]=lsim(G,u,t);
n=3;
D=zeros(N-n,n);
for i=n:N-1
    for j=0:n-1
        D(i-n+1,j+1)=y(i-j);
    end
end
Y=y(n+1:N,1);
a=((D'*D)^(-1))*D'*Y

syms x
eq=[x^3+a(1)*x^2+a(2)*x+a(3)==0];
X=vpa(solve(eq,x),2)
z=zeros(4,1);
for i=1:3
    z(i)=X(i);
end
Z=zeros(N,n);
for i=0:N-1
    Z(i+1,:)=[z(1)^i z(2)^i z(3)^i];
end
B=(Z'*Z)^(-1)*Z'*y
h=zeros(N,1);
for i=1:N
    h(i)=B(1)*z(1)^t(i)+B(2)*z(2)^t(i)+B(3)*z(3)^t(i);
end
%% Plot

plot(t,y,t,h)
grid on
legend('impulse response of system','impulse response of model')
xlabel('time(second)')
ylabel('h(t)')