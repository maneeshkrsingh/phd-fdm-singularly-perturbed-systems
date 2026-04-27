% BT-CS scheme for Heat equation
clear all
clc
L=input('Enter total length of spatial var')
T=input('Enter the upper bound for time')
h=input('Enter length of each spatial grid')
k=input('Enter length of each time grid')
alpha=input('Enter the duffison constant')
M=L/h;
N=T/k;
r=alpha*(k/h^2);
a=1+2*r;
b=-r;
c=-r;
A=zeros(M-1);
U=zeros(M+1,1);
D=zeros(M-1,1);
% Write boundry conditions in discrete forms
 
% U(1,1)=0;
%  U(M+1,1)=1;


 
for j=1:M-1
 
        D(j,1)=(((j-1)*h)^2);

end

for i=1:N+1
U(1,1)=(i-1)*k;
 U(M+1,1)=((((i-1)*k)^2)+1);    
D(1,1)=D(1,1)+(-c)*(i-1)*k;
D(M-1,1)=D(M-1,1)+(-b)*((((i-1)*k)^2)+1);

% Construct finite diffence scheme using matrix method

for k=1:M-1
    for j=1:M-1
        if(j==k)
            A(k,j)=a;
            if(j~=M-1)
            A(k,j+1)=b;
            A(k+1,j)=c;
            end
        end
    end
end
% A
% D(1,1)=k;
% D(M-1,1)=(((M-1)*k)^2)+1;
% D
% for i=1:N
%     for j=1:M
%     U(i+1,j)=inv(A)*(U(i,j))'+D;
%     U(i+1,1)=i*k;
%     U(i+1,M+1)=((i*k)^2)+1;
%     end
% end

U(2:M,1)=inv(A)*D;
D=U(2:M,1);
end
U
% surfc(0:h:L,0:k:T,U)
         
        
