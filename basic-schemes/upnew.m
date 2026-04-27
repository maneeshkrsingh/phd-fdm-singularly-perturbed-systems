clear all
clc
format long 
L=input('Enter total spatial size')
%h=input('Enter length of each spatial grid')
%e=input('Enter the value of epsilon')
for ep_counter=1:18
     e=2^(-ep_counter);
     for N_t_counter= 1:7
          N=2^(2+N_t_counter);
h=1/N;
r=2*h/e;
A=zeros(N-1);
U=zeros(N+1,1);
D=zeros(N-1,1);
a=-(2+r);
b=1+r;
c=1;
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
A;
D(1,1)=-1;
D(N-1,1)=0;
D;

% Construct upwind diffrence operator

 U(2:N,1)=inv(A)*D;
 U(1)=1;
 U(N+1)=0;
U;
%plot(0:h:1,U)
%hold on
y=0:h:1;
 U_E=(exp(-2*y/e) - exp(-2/e))/(1-exp(-2/e));
 
 err=abs(U-(U_E)');

err_sup=max(err);
%max_err=max(err_sup);
error(ep_counter,N_t_counter)=err_sup;
     end
end
err=abs(U-(U_E)');

err_sup=max(err);
%max_err=max(err_sup);
error(ep_counter,N_t_counter)=err_sup
    max_error=max(error)

%surf(ep_counter,N_t_counter,error)
%counter;

%U_E
%plot(y,U_E,'r--');
%hold off

