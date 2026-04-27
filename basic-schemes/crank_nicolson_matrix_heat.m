clear all
clc
L=input('Enter the total spatial size\n')
T=input('Enter the total time size\n')
del_t=input('Enter the size of each time grid\n')
del_x=input('Enter the size of each spatial grid\n')
alpha=input('Enter the value of diffusion constant\n')
N=L/del_x;
M=T/del_t;

r=alpha*del_t/(del_x*del_x);

A=zeros(N-1);
B=zeros(N-1);
U=zeros(M,N);
D=zeros(N-1,1);

    for i=1:N+1
         U(1,i)=((i-1)*del_x)^2;
    end

a=2-2*r;
b=r;
c=r;


for j=1:M
    for i=1:N-1
        if(i==j)
            A(j,i)=a;
            if(i~=N-1)
            A(j,i+1)=b;
            A(j+1,i)=c;
            end
        end
    end
end
A
d=2+2*r;
e=-r;
f=-r;

for j=1:M
    for i=1:N-1
        if(i==j)
            B(j,i)=d;
            if(i~=N-1)
            B(j,i+1)=e;
            B(j+1,i)=f;
            end
        end
    end
end
B
D(1,1)=del_t;
D(N-1,1)=(((N-1)*del_t)^2)+1;
D

for j=1:M
    U(j+1,2:N)=inv(B)*A*(U(j,2:N))'+D;
    U(j+1,1)=(j+1)*del_t;
    U(j+1,N+1)=((j+1*del_t)^2)+1;
end
U    
surf(0:del_x:L,0:del_t:T,U);
