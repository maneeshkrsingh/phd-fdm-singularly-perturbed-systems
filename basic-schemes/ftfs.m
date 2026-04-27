% FT-FS scheme for hyperbolic pde
clear all
clc
L=input('Enter the total length of spatial var')
T=input('Enter the upper bound for time grid')
h=input('Enter the length of each spatial grid')
k=input('Enter the length of each time grid')
M=L/h;
N=T/k;
R=(-k)/h;
% Enter the given Boundry conditions
for j=1:M
    U(1,j)=(sin(pi*(j*h)))^40;
end
for n=1:N
    U(n,1)=0;
    U(n,M+1)=0;
end
% Construct finite difference FT-FS scheme
for n=1:N-1
    for j=2:M
        U(n+1,j)=(1+R)*U(n,j)-R*U(n,j+1);
    end
end
for n=1:N-1
     U(n+1,M)=(1+R)*U(n,j);
end
for n=1:N-1
    U(n+1,2)=-R*U(n,2);
end
U
plot(0:h:L,U)
xlabel('spatial variable')
ylabel('value of u at spatial grids')

