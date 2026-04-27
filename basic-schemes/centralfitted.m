clear all
clc
L=input('Enter total spatial size')
h=input('Enter length of each spatial grid')
e=input('Enter the value of epsilon')
r=h/e;
N=L/h;
A=zeros(N-1);
U=zeros(N+1,1);
D=zeros(N-1,1);
a=-2;
b=1+r
c=1-r
U(1)=1;
U(N+1)=0;
for i=1:N
    for j=1:N-1
        if(j==i)
            A(i,j)=-2;
            if(j~=N-1)
                A(i,j+1)=1+r;
                A(i+1,j)=1-r;
            end
        end
    end
end
A
D(1,1)=-(1-r);
D(N-1,1)=0
D

% Construct central diffrence operator

 U(2:N,1)=inv(A)*D;
 U(1)=1;
 U(N+1)=0;
U
plot(0:h:1,U)
hold on
x=linspace(0,1,6);
v=exp((-2*x/e))-exp(-2/e);
w=1-exp(-2/e);
u=v/w;
plot(x,u,'--rs')
hold off