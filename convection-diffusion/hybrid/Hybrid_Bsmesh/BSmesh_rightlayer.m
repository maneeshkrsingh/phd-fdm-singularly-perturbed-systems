clc
clear all
format long e
N=16;
x=zeros(N+1,1);
%t=zeros(N+1,1);
h=zeros(N-1,1);
p=zeros(N-1,1);
e2=2^(-5);
e1=2^(-10);
%e3=2^(-8);
%sig3=2*e3*log(N);
tau=2.2;
sig2=tau*e2*log(N);
sig1=tau*(e1)*log(N);
x(1)=0;
% for i=1:(N/4)
%     x(i+1)=-2*e1*log(1-4*(1-(1/N))*(i/N));
% end
% for i=(N/4)+1:(N/2)
%     x(i+1)=x((N/4)+1)-2*(e2-e1)*log(1-4*(1-(1/N))*(i/N-1/4));
%     %x(i+1)=x((N/4)+1)+i*(4*((sig2)-(sig1)))/N;
% end

for i=1:(N/2)
    x(i+1)=(1-tau*e2*log(N))*2*i/N;
end
for i=(N/2)+1:(3*N/4)
    x(i+1)=(1-sig1)+tau*(e2-e1)*log(1-4*(1-(1/N))*(3/4-i/N));
    %x(i+1)=x((N/4)+1)+i*(4*((sig2)-(sig1)))/N;
end

for i=(3*N/4)+1:(N)
   x(i+1)=1+tau*e1*log(1-4*(1-(1/N))*(1-i/N));
end

% for i=(N/2)+1:N-1
%    x(i+1)=1-((1-2*e3*log(N))*2*(N-i)/N);
% end
x(N+1)=1;

x;
for i=1:N
h(i)=x(i+1)-x(i);
end
h;
for i=1:N-1
p(i)=h(i+1)-h(i);
end
x
h
plot(x)