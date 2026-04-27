
clear all
clc
L=input('Enter the total spatial size\n');
del_x=input('Enter the size of each spatial grid\n');
N=L/del_x;
e=input('Enter the value of epsilon');
a=2; %not fixed..will change depending on the problem

p=a/(N*e);
sig=(p/2)*coth(p/2);

A=zeros(N-1);
U=zeros(N+1,1);
B=zeros(N-1,1);

for j=1:N-1
    for i=1:N-1
        if(i==j)
             if(i==1)          %only giving for i is sufficient,bcoz i==j is given
                A(j,i)=-e*(sig+del_x*a/(2*e))/del_x^2;
             else
                A(j,i)=-2*e*sig/(del_x^2);
             end
            if(i~=N-1)
            A(j,i+1)=e*(sig+del_x*a/(2*e))/del_x^2;
            A(j+1,i)=e*(sig-del_x*a/(2*e))/del_x^2;
            end
        end
    end
end


B(1,1)=-(e*(sig-del_x*a/(2*e))/del_x^2)*(2*del_x/e);
B(N-1,1)=0;

U(1,1)=1;
U(N+1,1)=0;
U(2:N,1)=inv(A)*B
hold on
plot(0:del_x:1,U,'-+b');

y=0:0.001:1;
 U_E=(exp(-2*y/e) - exp(-2/e));


U_E
plot(y,U_E,'r--');
hold off
