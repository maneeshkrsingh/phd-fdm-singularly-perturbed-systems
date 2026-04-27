clear all
clc
x=input('Enter initial point of  spatial variable\n')
y=input('Enter final point of spatial variable\n')
T=input('Enter the upper bound for time grid\n')
h=input('Enter length of each spatial grid\n')
k=input('Enter size of each time grid\n')
alpha=input('Enter the value of diffusion constant\n')
M=(y-x)/h;
N=T/k;
r=alpha*(k/h^2);
a=1-2*r;
b=r;
c=r;
M
N
A=zeros(M-1);
U=zeros(N,M+1);
B=zeros(M-1,1);
% Construct Boundry Conditions
for j=1:N
    U(j,1)=(j*k)^6;
    U(j,M+1)=(j*k)^3+1;
end
for i=2:M
    U(1,i)=(i*h)^2
end
% Construct matrix A and B
for j=1:N
    for i=1:M-1
        if(i==j)
            A(j,i)=a;
            if(i~=M-1)
                A(j,i+1)=b;
                A(j+1,i)=c;
            end
        end
    end
end
B(1,1)=r*k^6;
B(M-1,1)=r*((M-1)*k)^3+1;

for j=1:N
    U(j+1,2:M)=A*(U(j,2:M))'+B;
    U(j+1,1)=0;
    U(j+1,M+1)=0;
    
end
A;
U    
surfc(0:h:y-x,0:k:T,U)
    
