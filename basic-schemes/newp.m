clear all
clc
L=input('Enter total spatial size')
for l=1:18
    for m=1:7
        e(1,l)=2^(-l);
        N(m,1)=2^(2+m);
    end
end
e
N

h=L./N;
r=2.*h./e;
A=zeros(N-1);
U=zeros(N+1,1);
D=zeros(N-1,1);
a=-(2+r);
b=1+r
c=1
U(1)=1;
U(N+1)=0;
for i=1:N
    for j=1:N-1
        if(j==i)
            A(i,j)=-(2+r);
            if(j~=N-1)
                A(i,j+1)=1+r;
                A(i+1,j)=1;
            end
        end
    end
end
D(1,1)=-1;
D(N-1,1)=0
% Construct central diffrence operator

 U(2:N,1)=inv(A)*D;
 U(1)=1;
 U(N+1)=0;
U
plot(0:h:1,U)
hold on
x=linspace(0,1,N);
v=exp((-2*x/e))-exp(-2/e);
w=1-exp(-2/e);
u=v/w;
error=abs(u'-U)
