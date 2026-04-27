% FT-BS scheme for Hyperbolic pde
N=input('Enter a number')
M=input('Enter a number of interval along time')
T=input('Enter the upper bound for time')
h=1/N;
k=T/M;
% Construct mesh point in spatial variable
x(1)=0;
for i=1:N-1
    x(i+1)=x(i)+h;
end
% construct mesh point in time variable
t(1)=0;
for j=1:M
    t(j+1)=t(j)+k;
end
% Generate time leval atj=1
for i=1:N
    U(i,1)=2*x(i)+1;
end
% Generate function at i=1
for j=1:M
   U(1,j+1)=5;
end
% Define the CFL quantity
R=k/h
for j=1:M-1
    for i=2:N
        U(i,j+1)=(1-R)*U(i,j)+R*U(i-1,j);
    end
end

U
surfc(0:h:1,0:k:T,U)