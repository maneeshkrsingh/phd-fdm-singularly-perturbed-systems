
clear all
clc
% FT-CS scheme for Heat equation
h=input('Enter size of each spatial grid\n')
k=input('Enter size of each time grid\n')
T=input('Enter upper bound for time interval\n')
M=1/h;
N=T/k;
r=(k/h^2);
%construct grid points along spatial axis
for j=1:M
    x(j+1)=j*h;
end
for i=1:N
    t(i)=i*k;
end
% write boundry conditions in discrete form
for j=1:M+1
    U(1,j)=(j*h)^2;
end
for i=1:N
    U(i,1)=0;
end
for i=1:N
    U(i,M+1)=(((i)*k)^5)+1;
end
% Construst FT-CS difference schene
for i=1:N
    U(i+1,2)=(1-2*r)*U(i,2)+r*U(i,3);
end
for i=1:N
    for j=3:M-1
        U(i+1,j)=(1-2*r)*U(i,j)+r*U(i,j+1)+r*U(i,j-1);
    end
end
for i=1:N
    U(i+1,M)=(1-2*r)*U(i,M)+r*U(i,M-1)+(((i+1)*k)^5)+1;
end
U
surfc(0:h:1,0:k:T,U)
   