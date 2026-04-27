function [k,x]=polation(alfasolve,nel)
syms k
%alfasolve=rand(3*(N+1),1);
x(1)=0;
h=1/nel;
N=nel;

for i=2:N+1
    x(i)=x(1)+(i-1)*h;
end
y = sym('y');
% phi_0=@(y)1;
% phi_1=@(y)2*(y-(x(i)+x(i+1))/2)/(x(i+1)-x(i));
% phi_2=@(y)4*(y-(x(i)+x(i+1))/2)/(x(i+1)-x(i))^2;

i=0;
sum=0;
for j=1:3:3*N
    i=i+1;
    if(i<=N)
%     phi(j)=1;
%     phi(j+1)=2*(y-(x(i)+x(i+1))/2)/(x(i+1)-x(i));
%     phi(j+2)=4*(y-(x(i)+x(i+1))/2)^2/(x(i+1)-x(i))^2;
sum=sum+1*alfasolve(j)+2*(y-(x(i)+x(i+1))/2)/(x(i+1)-x(i))*alfasolve(j+1)+4*(y-(x(i)+x(i+1))/2)^2/(x(i+1)-x(i))^2*alfasolve(j+2);
%sum1=1*alfasolve(j)+2*(y-(x(i)+x(i+1))/2)/(x(i+1)-x(i))*alfasolve(j+1)+4*(y-(x(i)+x(i+1))/2)^2/(x(i+1)-x(i))^2*alfasolve(j+2)
    else
        break;
    end
end
k=symfun(sum,y);

