% BT-FS scheme for hyperbolic pde
clear all
clc
L=input('Enter the total lenth of spatial variable')
T=input('Enter the upper bound for time grid')
h=input('Enter the length of each spatial grids')
k=input('Enter the length of each time grids')
N=L/h;
M=T/k;
R=k/h;

A=zeros(N);
U=zeros(N+1,1);
B=zeros(N,1);
% write given boundry conditions
for j=1:N
    B(j,1)=1;
end
U(N+1,1)=1; 
% Costruct matrix A 
for i=1:M
    for j=1:N
        if (j==i)
            A(i,j)=1-R;
             if(j<=N-1)
                A(i,j+1)=R;
             end
        end
    end
end
A

U(1:N)=inv(A)*B
plot(0:h:1,U)