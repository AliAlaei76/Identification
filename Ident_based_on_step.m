%% Continuous state space representation and transfer function of system

A=[0 1 0;-77.3707 -0.2703 -6.2447;-0.9879 0 -0.322];
B=[0 0 1].';
C=[1 0 0];
D=0;
sys=ss(A,B,C,D); 
x0=[0.5468 0 9.6926].';
s=tf('s');
G=zpk(sys)
%% step response (zero and nonzero initial condition)

t=0:0.01:40;
N=numel(t);
u=ones(1,N);
[y1,t]=lsim(G,u,t);
[y2,t]=lsim(sys,u,t,x0);
plot(t,y1,t,y2)
grid on
legend('Step response in zero initial condition','Step response in nonzero initial condition')
plot(t,y1)
grid on
legend('Step response in zero initial condition')
plot(t,y2)
grid on
legend('Step response in nonzero initial condition')
%% System identification based on step response of zero initial condition

k=-0.3331;
tau=4.01;
G2=k/(1+tau*s);
[y3,t]=lsim(G2,u,t);
plot(t,y1,t,y3,'linewidth',1.5);
grid on
legend('step response of system','step response of model')
xlabel('Time(second)')
ylabel('S(t)')
%% System identification based on step response of nonzero initial condition

M=-1.993/k-1
T=1.07-0.36
zeta=(-log(M))/(pi^2+(log(M)^2))^0.5
wn=(2/T)*(pi^2+(log(M)^2))^0.5
zeta=0.0199;
wn=8.796;
G3=(k*wn^2)/(s^2+2*zeta*wn*s+wn^2)
step(G3)