clear all
clc
x=linspace(0,1,6)
e=input('Enter the value of epsilon')
v=exp((-2*x/e))-exp(-2/e);
w=1-exp(-2/e);
u=v/w;
plot(x,u)