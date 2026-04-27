clc
clear all
x1=2*rand(1,100)-1;
y1=2*rand(1,100)-1;
z1=3./(1+x1.^2+y1.^2);
x2=linspace(-1,1,30);
y2=x2';
[Xi,Yi,Zi]=griddata(x1,y1,z1,x2,y2,'nearest');
