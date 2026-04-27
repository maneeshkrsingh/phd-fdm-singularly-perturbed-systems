clear all
clc
e=input('Enter the value of epsilon')
x=linspace(0,0.001);
v=50*(1-exp(-x/e));
w=1-exp(-1/e);
u=v/w;
plot(x,u,'b--')
hold on
plot(x,50,'b')
hold off