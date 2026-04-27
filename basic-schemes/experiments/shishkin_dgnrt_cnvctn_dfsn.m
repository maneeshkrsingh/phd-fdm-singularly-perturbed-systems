clc
clear all;
format short e
%  N=input('enter the value of N=');
%  N0=input('enter the value of N0=');
% e=input('enter the value of e=');

p=5;
for l=7:7

e(l)=1/(2^(5*(l-1)));
% e(l)=1/(10^(l-1));
for l2=1:2
    N=16*2^(l2-1);
    N0=16*2^(l2-1);
    
s=min(1/2,(e(l)^(1/2))*log(N));

%%%%%%%%%%%%%%%%%%
    N1=N;
    N01=N0;
%%%%%%%%%%%%%%%%%
% temp=1024/N;
for l1=1:2
    N=N1*2^(l1-1);
    N0=N01*2^(l1-1);
% if(l1>1)
%     N=1024;
%     N0=1024;
% end

A=zeros(N-1);
b=zeros(N-1,1);
d1=zeros(N-1,1);
u=zeros(N+1,N0+1);
 x=zeros(N+1,1);
 h=zeros(N,1);
t=zeros(N0+1,1);

% s=min(1/4,e(l)^(1/2)*log(N));

for i=1:(N/2)+1
    x(i)=2*s*(i-1)/N;
   
end

for i=(N/2)+2:N+1
    x(i)=s+(2*(i-1-(N/2))*(1-s)/N);
    
end
x;
for i=1:N
    h(i)=x(i+1)-x(i);
end
h;
%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:N0+1
    t(j)=(1*(j-1))/N0;
end
k=1/N0;
%%%%%%%%%%%%%%%%%
t;
%%%%%%%%%%%%%%%%%%%%%%%%%
% t(1)=0;
% k=0.8/N0;
% for j=1:N0
%     t(j+1)=t(j)+k;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N+1
    u(i,1)=(1-x(i))^2;
%       u(i,1)=(1-x(i))^2+1;
%       u(i,1)=0;
end
for j=1:N0+1
    u(1,j)=1+(t(j)^2);
%       u(1,j)=2+(t(j)^2);
%       u(1,j)=0;
    u(N+1,j)=0;
%     u(N+1,j)=1;
%      u(Nx+1,j)=-(t(j)^3)/3;
end


 %check it Nt+1
 for j=2:N0+1 
for i=2:N-2
        
   
       A(i,i-1)=(2*e(l))/((h(i)+h(i+1))*h(i));
       A(i,i)=-(((x(i+1))^p)/h(i+1))-(2*e(l))/((h(i)+h(i+1))*h(i+1))-(2*e(l))/((h(i)+h(i+1))*h(i))-(1/k)-1;
       A(i,i+1)=(2*e(l))/((h(i)+h(i+1))*h(i+1))+(((x(i+1))^p)/h(i+1));
   
       
end

i=1;
      A(i,i)=-(((x(i+1))^p)/h(i+1))-(2*e(l))/((h(i)+h(i+1))*h(i+1))-(2*e(l))/((h(i)+h(i+1))*h(i))-(1/k)-1;
       A(i,i+1)=(2*e(l))/((h(i)+h(i+1))*h(i+1))+(((x(i+1))^p)/h(i+1));
i=N-1;     
     A(i,i-1)=(2*e(l))/((h(i)+h(i+1))*h(i));
       A(i,i)=-(((x(i+1))^p)/h(i+1))-(2*e(l))/((h(i)+h(i+1))*h(i+1))-(2*e(l))/((h(i)+h(i+1))*h(i))-(1/k)-1;
    

   A;
   
% for j=2:N0+1   %check itNt+1
         for i=2:N-2
%              u(i+1,j-1)
              
                   b(i)=(x(i+1)^2)-1-((1/k)*u(i+1,j-1));
%                  b(i)=-(2*t(j))+(t(j)^2)-x(i+1)/2;
%                      b(i)=(((x(i+1)^2)+1)*(t(j)^2))-((1/k)*u(i+1,j-1));
              
         end       

         b(1)=(x(2)^2)-1-((1/k)*u(2,j-1))-((2*e(l))/((h(1)+h(2))*h(1)))*u(1,j);
          b(N-1)=(x(N)^2)-1-((1/k)*u(N,j-1))-((2*e(l))/((h(N-1)+h(N))*h(N))+(((x(N))^p)/h(N)))*u(N+1,j);
%          b(N-1)=-(2*t(j))+(t(j)^2)-x(N)/2;
%           b(1)=(((x(2)^2)+1)*(t(j)^2))-((1/k)*u(2,j-1))-((2*e(l))/((h(1)+h(2))*h(1)))*u(1,j);
%           b(N-1)=(((x(N)^2)+1)*(t(j)^2))-((1/k)*u(N,j-1))-((2*e(l))/((h(N-1)+h(N))*h(N))+(((x(N))^p)/h(N)))*u(N+1,j);
         
          d1=A\b;
%           d=inv(A)*b;

         for i=2:N
              u(i,j)=d1(i-1); %storing the value of d(each time level) in column in u.
         end

 end
      e(l);
      if(l1==1)
          F=u;
      end
       if(l1==2)
          G=u;
      end
end

%      N
%      N0
%      max_err(l,l2)=max(max(abs(F-G(1:temp:N+1,1:temp:N0+1))));
        max_err(l,l2)=max(max(abs(F-G(1:2:N+1,1:2:N0+1))));
        max_err
         if(l2>1)
         cnv_rt(l,l2-1)=log2(max_err(l,l2-1)/max_err(l,l2))
         end
end
end
% A
% b
% surf(x,t,u')
% max_max_err=max(max_err)
% for i=1:5
%     cnv_rt(i)=log2(max_max_err(i)/max_max_err(i+1))
% end