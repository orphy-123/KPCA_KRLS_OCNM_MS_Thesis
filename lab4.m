clc
clear all
t=0:0.01:1;
f=1;
x=2*sin(2*pi*f*t);

fs=10; % Sampling Frequency
Ts=1/fs;
t1=0:Ts:1;
xs=2*sin(2*pi*f*t1);
figure(1)
plot(t,x);
hold on
plot(t1,xs,'r');
figure(2)
stem(t1,xs,'r')


partition=[-2:0.5:2]; %9 levels
codebook=[-2,-1.75:0.5:1.75,2] %11 assigned points
[index,xq]=quantiz(xs, partition, codebook);
figure(3)
plot(t1,xs,'x',t1,xq,'.')
legend('Sampled Signal','Quantized Signal')

e=abs(xs-xq);
l=length(e);
s=sum(e.^2);
q_noise_power=s/l