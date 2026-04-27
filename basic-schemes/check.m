clear all
clc
format long 
L=input('Enter total spatial size')
sigo=input('Enter the value of sigma zero')

for ep_counter=1:10
     e=2^(-2*(5+ep_counter));
    % for N_t_counter= 1:6
          N=2^(6);
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
% a(1)=func(0);
 for j=1:N
  a(1)=1;   
a(j)=func(x(j));
a(N+1)=1;
 end
 % construct hodie coeffcient
 for j=1:N/2
     q(j)=a(j+1)/(a(j+1)+a(j));
 end
 for j=(N/2)+1:N-1
     q(j)=(h(j)+h(j+1))/(3*h(j));
 end
 % construct finite dfference scheme
A=zeros(N-1);
U=zeros(N+1,1);
D=zeros(N-1,1);
for j=1:N-1
    l(j)=((-2*e)+q(j)*(-((2*h(j)+h(j+1))*a(j)))-(1-q(j))*h(j+1)*a(j+1))/(h(j+1)*(h(j)+h(j+1)));
    n(j)=((-2*e)+(h(j)*a(j+1))-(q(j)*h(j)*a(j+1)+a(j)))/(h(j+1)*(h(j)+h(j+1)));
    m(j)=-l(j)-n(j);
    
    for k=1:N-1
        if(k==j)
            A(j,k)=m(j);
            if(k~=N-1)
                A(j,k+1)=n(j);
                A(j+1,k)=l(j);
            end
        end
    end
    A;
    inv(A);
D(1,1)=-l(1)+nonh(x(2),e);
D(N-1,1)=nonh(x(N),e);
for i=2:N-2
    D(i,1)=nonh(x(i+1),e);
end
D;

U(2:N,1)=inv(A)*D;
 U(1)=1;
 U(N+1)=0;
%  plot(x,U)

end
x;
U;
U_E=(1-exp(-(1-x)/e))/(1-exp(-1/e))-cos(pi*x/2);
err=abs(U-(U_E)')

err_sup=max(err);
max_err=max(err_sup);
error(ep_counter,N)=max_err;

    
 
     %end
end
 max_error=max(error)
       