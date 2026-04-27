clc
clear all
format long

for ep_counter=1:10
    e=2^(-(ep_counter-1));
 for N_counter=1:4
     N=4^(N_counter+1);
     
     sig=min(1/4,sqrt(e)*log(N));
     
for i=1:(N/4)
    h(i)=(4*sig)/N;
end
for i=(N/4)+1:3*N/4
    h(i)=(2*(1-2*sig))/N;
 end
for i=(3*N/4)+1:N
    h(i)=(4*sig)/N;
end
 for i=1:N
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
x;  
 for l=1:2 
     if (l==2)
         N=2*N;
         sig=min(1/4,sqrt(e)*log(N/2));
for i=1:(N/4)
    h(i)=(4*sig)/N;
end
for i=(N/4)+1:3*N/4
    h(i)=(2*(1-2*sig))/N;
 end
for i=(3*N/4)+1:N
    h(i)=(4*sig)/N;
end
 for i=1:N
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
x;  
     end
         
     

A=zeros(N-1);
U=zeros(N+1,1);
D=zeros(N-1,1);


 
  U(1)=0;
  U(N+1)=0;
for i=1:N+1
    f(i)=x(i)*(1-x(i));
end   

 for i=2:N-2
     A(i,i-1)=-2*e/((h(i)+h(i+1))*h(i));
     A(i,i)=2*e/(h(i)*h(i+1))+1+x(i);
     A(i,i+1)=-2*e/((h(i)+h(i+1))*h(i+1));
 end
 i=1;
       A(i,i)=2*e/(h(i)*h(i+1))+1+x(i);
     A(i,i+1)=-2*e/((h(i)+h(i+1))*h(i+1));
     
 i=N-1;   
  A(i,i-1)=-2*e/((h(i)+h(i+1))*h(i));
     A(i,i)=2*e/(h(i)*h(i+1))+1+x(i);

 for i=1:N-1
     D(i)=f(i+1);
 end
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

error
cnv
% U(2:N,1)=inv(A)*D';
%   plot(x,U);

