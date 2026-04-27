clc
clear all
format long e
N=32;
x=zeros(N+1,1);
%t=zeros(N+1,1);
h=zeros(N-1,1);
p=zeros(N-1,1);

 e1=2^(-10);
    e2=2^(-6);
    sig2=2*e2*log(N);
sig1=2*(e1)*log(N);
x(1)=0;
for i=1:(N/4)
    x(i+1)=2*e1*((i/N)*log(N))/(0.25+(0.25-(i/N))*log(N));
end
for i=(N/4)+1:(N/2)
    x(i+1)=x((N/4)+1)+2*(e2-e1)*((((i/N)-0.25)*log(N))/(0.25+(0.25-((i/N)-0.25))*log(N)));
    %x(i+1)=x((N/4)+1)+i*(4*((sig2)-(sig1)))/N;
end
% for i=(N/2)+1:N-1
%    x(i+1)=2*(1-(2*e2*log(N)))*(i-1)/N;
% end

for i=(N/2)+1:N-1
   x(i+1)=1-((1-2*e2*log(N))*2*(N-i)/N);
end
x(N+1)=1;

x;
for i=1:N
h(i)=x(i+1)-x(i);
end
h;
for i=1:N-1
p(i)=h(i+1)-h(i);
end

%plot(h,N)