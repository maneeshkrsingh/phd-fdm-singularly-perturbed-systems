clear all
clc
% Euler Scheme
N=input('Enter a number')
M=input('Enter a number of intervals along time')
T=input('Enter the upper bound for time')
h=1/N;
x(1)=0;
for i=1:N
    x(i+1)=x(i)+h;
end
k=T/M
t(1)=0;
for j=1:M
    t(j+1)=t(j)+k;
end
% Generate the known time leval j=1
for j=1:M
    for i=1:N+1
    if(j==1)
        if((i-1)*h<=0.5)
    U(i,j)=2*(i-1)*h;
        end
        if((i-1)*h>=0.5)
            U(i,j)=2*(1-(1-i)*h);
        end
    end
    end
end
% define CFL quantity where stability depends
r=k/h^2;
for j=1:M
    for i=2:N
        U(i,J+1)=r*U(i-1,j)+(1-2*r)*U(i,j)+r*U(i+1,j);
    end
end
U
surfc(0:h:1,0:k:T,U)