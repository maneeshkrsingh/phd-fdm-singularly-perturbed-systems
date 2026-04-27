clear all
clc
x=linspace(0,0.001,10);
e=input('enter the value of epsilon')
v=50*(1-exp(-x/e));
w=(1-exp(-1/e));
u=v/w;
plot(x,u)