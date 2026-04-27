clear all
clc
L=input('Enter the total length of spatial var')
T=input('Enter the upper bound of time')
h=input('Enter the length of each special grid')
k=input('Enter the length of each time grid')
alpha=input('Enter the diffusion constant')
N=L/h;
M=T/k;
r=alpha*(k/h^2);
a=1+2*r
b=-r
c=-r
A=zeros(N-1);
U=zeros(M,N);
D=zeros(N-1,1);
% Write initial and boundry conditions in discrete forms
for i=1:N+1
    U(1,i)=((i-1)*h)^2;
end
for j=1:M
    U(j,1)=(j-1)*k;
    U(j+1,N+1)=(((j-1)*k)^2)+1;
end
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
D(1,1)=k;
D(N-1,1)=(((N-1)*k)^2)+1;
D

for j=1:M
    U(j+1,2:N)=inv(A)*(U(j,2:N))'+D;
    U(j+1,1)=(j)*k;
    U(j+1,N+1)=((j*k)^2)+1;
end
U
surf(0:h:L,0:k:T,U)