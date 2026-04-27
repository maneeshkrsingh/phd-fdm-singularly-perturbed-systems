clc
clear all
format short e
N=64;
x=zeros(N+1,1);
%t=zeros(N+1,1);
h=zeros(N-1,1);
p=zeros(N-1,1);

 e1=2^(-10);
 e2=2^(-6);
 
sig2=1.2*e2*log(N);
sig1=1.2*e1*log(N);
x(1)=0;


for i=1:(N/4)
    x(i+1)=4*(0.5-sig2)*i/N;
end
for i=(N/4)+1:(3*N/8)
    %x(i+1)=x((N/4)+1)+1.2*(e2-e1)*log(1-8*(1-1/N)*(3/8-i/N));
     x(i+1)=(0.5-sig1)+1.2*(e2-e1)*log(1-8*(1-1/N)*(3/8-i/N));
   
end

for i=(3*N/8)+1:(N/2)
   %x(i+1)=x((3*N/8)+1)+1.2*e1*log(1-8*(1-1/N)*(4/8-i/N));
    x(i+1)=0.5+1.2*e1*log(1-8*(1-1/N)*(4/8-i/N));
end

for i=(N/2)+1:5*N/8
   x(i+1)=0.5-1.2*e1*log(1-8*(1-1/N)*(i/N-1/2));
end


for i=(5*N/8)+1:6*N/8
   %x(i+1)=x((5*N/8)+1)-1.2*(e2-e1)*log(1-8*(1-1/N)*(i/N-5/8));
   x(i+1)=(0.5+sig1)-1.2*(e2-e1)*log(1-8*(1-1/N)*(i/N-5/8));
end

for i=(6*N/8)+1:N-1
   x(i+1)=(0.5+sig2)+(i-6*N/8)*(0.5-sig2)*4/N;%1-((1-2*e3*log(N))*2*(N-i)/N);
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