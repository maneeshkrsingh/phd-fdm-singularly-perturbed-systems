clc;
clear all;
close all;
HTOL=0.001;
dH=1;
N=2;
% h=1/N;
%e=input('Enter the value of singulaer perturbation parameter  ');
 e=0.001;
 itrn=0;
%   for m=1:4
   
A=zeros(N-1);
b=zeros(N-1,1);
B=zeros(N-1);
f=zeros(N-1,1);
M=zeros(N,1);
u(1)=0;
u(N+1)=1;

for i=1:N+1
x(i)=(i-1)/N;
end

while(dH>HTOL )
  itrn=itrn+1;
for i=1:N
h(i)=x(i+1)-x(i);
end

   for i=1:N-1
    for j=1:N-1
        if(i==j)
            %A(i,j)=((2*e)/(h(j)*h(j+1)))+1/h(j);
            A(i,j)=((2*e)/(h(i)*h(i+1)))+1/h(i);
        end
            if(j==i-1)
%                 A(i,j)=-((2*e)/(h(j+1)*(h(j+2)+h(j+1))))-1/h(j+1);
                  A(i,j)=-((2*e)/(h(i)*(h(i)+h(i+1))))-1/h(i);
            end
                if(j==i+1)
%                  A(i,j)=-(2*e)/(h(j)*(h(j)+h(j-1))); 
                   A(i,j)=-(2*e)/(h(i+1)*(h(i)+h(i+1))); 
                end
    end
   end
A;
b(1,1)=0;
b(N-1,1)=(((2*e)/(h(N)*(h(N-1)+h(N)))))*u(N+1);
b;
 u(2:N)=A\b;
% u(2:N)=linsolve(A,b);

alpha=1;
for j=2:N+1
M(j-1)=sqrt(alpha + ((u(j)-u(j-1))/(x(j)-x(j-1)))^2);
H(j-1)=M(j-1)*(x(j)-x(j-1));
end
M;
dH=log10(max(H)/min(H));

    for i=1:N-1
    for j=1:N-1
        if(i==j)
%             B(i,j)=-M(j)-M(j+1);
              B(i,j)=-M(i)-M(i+1);
        end
            if(j==i-1)
%                 B(i,j)=M(j+1);
                  B(i,j)=M(i);
            end
                if(j==i+1)
%                  B(i,j)=M(j); 
                   B(i,j)=M(i+1); 
                end
    end
    end
B;
f(N-1,1)=-M(N)*x(N+1);
% x(1)=0;
% x(N+1)=1;
 y(2:N)=B\f;
% y(2:N)=linsolve(B,f);
y(1)=0;
y(N+1)=1;
x=y;

plot(x,itrn,'b-*')
hold on;
end
hold off;
x;
 u;
 for i=1:N+1
% u_exact(i)=(exp((x(i)-1)/e)-exp((-1)/e))/(1-exp(-1/e));
 u_exact(i)=(exp((x(i)-1)/e)-exp(-1/e))/(1-exp(-1/e));
 end
%  u_exact;
 for i=1:N+1
       err(i)=abs(u_exact(i)-u(i));
 end
%  err=abs(u_exact-u)
 err_max=max(err)
%  if(m>1)
%     conv_rate(m)=abs(log2(err_max(m)/err_max(m-1)));
%  end
%   
%  N=N*2;
%  dH=1;
%  HTOL=0.001;
%   end
%       err_max
%       conv_rate
% 
plot(x,u,'b--')
 hold on
plot(x,u_exact,'r--')
hold off
% legend('approximate','exact')
%  plot(x,err,'b-o')