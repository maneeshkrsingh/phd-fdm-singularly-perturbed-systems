clear all
clc
format short e
L=input('Enter total spatial size')
sigo=input('Enter the value of sigma zero')

for ep_counter=1:10
     e=2^(-2*(5+ep_counter));
     for N_t_counter= 1:6
          N=2^(5+N_t_counter);
          x=zeros(1,N+1);
 sig=sigo*e*log(N);
for i=1:N/2
   h(i)=2*(1-sig)/N;
end

for i=N/2+1:N
   h(i)=2*sig/N;
end

x(1)=0;
for i=1:N
   x(i+1)=x(i)+h(i);
end
a(1)=func(0);
%  for j=1:N+1
% a(j)=func(x(j))
%  end
     end
end
          