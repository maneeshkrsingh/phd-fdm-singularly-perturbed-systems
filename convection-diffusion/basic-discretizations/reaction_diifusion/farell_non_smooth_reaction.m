clc
clear all
format short

for ep_counter=1:5
    e=10^(-(ep_counter));
 for N_counter=1:2
     N=8^(N_counter);
  d=1/2;   
     sig1=min(d/4,2*sqrt(e)*log(N));
     sig2=min((1-d)/4,2*sqrt(e)*log(N));
     %sig1=min(1/8,sqrt(e)*log(N));
     %sig2=min(1/8,sqrt(e)*log(N));
     
for i=1:(N/8)
    h(i)=(8*sig1)/N;
end
for i=(N/8)+1:3*N/8
    h(i)=4*(d-2*sig1)/N;
 end
for i=(3*N/8)+1:N/2
    h(i)=(8*sig1)/N;
end

for i=N/2+1:5*(N/8)
    h(i)=(8*sig2)/N;
end
for i=(5*N/8)+1:7*N/8
    h(i)=4*((1-d)-2*sig2)/N;
 end
for i=(7*N/8)+1:N
    h(i)=(8*sig2)/N;
end

 for i=1:N
    x(1)=0;
    x(i+1)=x(i)+h(i);
 end
 for l=1:2 
     if (l==2)
         N=2*N;
         sig1=min(d/4,2*sqrt(e)*log(N/2));
     sig2=min((1-d)/4,2*sqrt(e)*log(N/2));
 
 for i=1:(N/8)
    h(i)=(8*sig1)/N;
end
for i=(N/8)+1:3*N/8
    h(i)=4*(d-2*sig1)/N;
 end
for i=(3*N/8)+1:N/2
    h(i)=(8*sig1)/N;
end

for i=N/2+1:5*(N/8)
    h(i)=(8*sig2)/N;
end
for i=(5*N/8)+1:7*N/8
    h(i)=4*((1-d)-2*sig2)/N;
 end
for i=(7*N/8)+1:N
    h(i)=(8*sig2)/N;
end

 for i=1:N
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
     end
A=zeros(N-1);
U=zeros(N+1,1);
D=zeros(N-1,1);

  U(1)=0;
  U(N+1)=0;
for i=1:N/2-1
     f(i)=0.7;
end 
 for i=N/2+1:N
     f(i)=-0.6;
end 
% for i=1:N/2
%      b(i)=(1+x(i)^2);
% end 
%  for i=N/2+1:N
%      b(i)=(2+x(i)^3);
% end

 for i=2:N/2-1
     A(i,i-1)=-2*e/((h(i)+h(i+1))*h(i));
     A(i,i)=2*e/(h(i)*h(i+1))+1;
     A(i,i+1)=-2*e/((h(i)+h(i+1))*h(i+1));
 end
 
 for i=N/2+1:N-2
     A(i,i-1)=-2*e/((h(i)+h(i+1))*h(i));
     A(i,i)=2*e/(h(i)*h(i+1))+1;
     A(i,i+1)=-2*e/((h(i)+h(i+1))*h(i+1));
 end
 
 i=1;
       A(i,i)=2*e/(h(i)*h(i+1))+1;
     A(i,i+1)=-2*e/((h(i)+h(i+1))*h(i+1));
 
  i=N/2;
     A(i,i-1)=1/h(i);
     A(i,i)=-(1/h(i)+1/h(i+1));
     A(i,i+1)=1/h(i+1);   
     
     
 i=N-1;   
     A(i,i-1)=-2*e/((h(i)+h(i+1))*h(i));
     A(i,i)=2*e/(h(i)*h(i+1))+1;

 for i=1:N/2-1
     D(i)=f(i);
 end
 for i=N/2+1:N-1
     D(i)=f(i);
 end
 i=N/2;
 D(i)=0;
%  if(i~=1 && i~=N-1)
%         if(i~=N/2)
%              D(i)=f(i+1);
%  %B(k,1)=-9;
%         end
%         if(i==N/2)
%             D(i)=0;
%         end
%      end
%      if(i==1)
%          D(i)=f(i+1);
%      end
%      if(i==N-1)
%         D(i)=f(i+1);
%      end
 U(2:N,1)=inv(A)*D;
if (l==1)
 U_pre=U;
else
    err=abs(U(1:2:N+1)-U_pre);
    max_err=max(err);
end

 end
  error(ep_counter,N_counter)=max_err;
 end

 end

for j=1:ep_counter
    for i=1:N_counter-1
         cnv(j,i)=log2(error(j,i)/error(j,i+1));
    end
end
plot(x,U)
 error
 cnv