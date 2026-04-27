clear all
clc
L=input('Enter the toltal length of spatial var')
T=input('Enter the upper bound of time')
h=input('Enter the length of each spatial grid')
k=input('Enter the length of each time grid')
a=input('Enter the speed of solution')
R=a*k/h;
N=L/h;
M=T/h;
A=zeros(N-1);
U=zeros(N,M);
D=zeros(N-1,1);
% Enter the given boundry conditions in discrete form
for i=1:N+1
    U(1,i)=((i-1)*h)^2;
end
for j=1:M
    U(j,1)=(j-1)*k;
    U(j,N+1)=(((j-1)*k)^2)+1;
end
% construct matrices required for writing diffence scheme
for i=1:N-1
    for j=1:N-1
        if(j==i)
            A(i,j)=(1+R^2);
            if(j~=N-1)
                A(i,j+1)=(((-R^2)/2)+R/2);
                A(i+1,j)=(((-R^2)/2)-R/2);
            end
        end
    end
end
A
D(1,1)=(((-R^2)/2)-R/2)*j*k;
D(N-1,1)=(((-R^2)/2)+R/2)*(((j*k)^2)+1);
D
% Construct laxwendroff scheme
for j=1:M
    U(j+1,2:N)=inv(A)*(U(j,2:N))'+D;
    U(j+1,1)=(j)*k;
    U(j+1,N+1)=((j*k)^2)+1;
end
U
plot(0:h:L,U)