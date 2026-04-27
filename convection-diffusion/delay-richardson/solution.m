
function [u_sol]= solution(e,N,p,x,del_t,u_d,time_counter)   
A=zeros(N-1);
B=zeros(N-1,1);
u=zeros(p+1,N+1);
u_sol=zeros(p+1,N+1);
U=zeros(N+1,1);
% D=zeros(N-1);
    u(1,:)=u_d(1,:);
    
for i=2:N+1
    j=i-1;
    h(j)=x(i)-x(i-1);
end
h;
for n=1:p
    x;
    for j=2:N
    k=j-1;
    f=(x(j+1)-x(j));
%     h(j)
    b=(x(j)-x(j-1));
%     h(j-1)
    c=(x(j+1)-x(j-1));
    if(k~=1 && k~=N-1)
        A(k,k)=(2*e*del_t/(f*c))+1+(2*e*del_t/(b*c))+(2-x(j)*x(j))*del_t/b+del_t*(x(j)+1)*((n*del_t+(time_counter-1))+1);
        A(k,k-1)=-(2*e*del_t/(b*c))-(2-x(j)*x(j))*del_t/b;
        A(k,k+1)=-(2*e*del_t/(f*c));
    end
    if(k==1)
        A(k,k)=(2*e*del_t/(f*c))+1+(2*e*del_t/(b*c))+(2-x(j)*x(j))*del_t/b+del_t*(x(j)+1)*((n*del_t+(time_counter-1))+1);
        A(k,k+1)=-(2*e*del_t/(f*c));
    end
    if(k==N-1)
        A(k,k)=(2*e*del_t/(f*c))+1+(2*e*del_t/(b*c))+(2-x(j)*x(j))*del_t/b+del_t*(x(j)+1)*((n*del_t+(time_counter-1))+1);
        A(k,k-1)=-(2*e*del_t/(b*c))-(2-x(j)*x(j))*del_t/b;
    end        
        
    end
m=n-1;
 n;
 p;
 del_t;
u(n+1,1)=0; %Boundary condition
u(n+1,N+1)=0; %Boundary condition
u;
B(1,1)=u(n,2)+del_t*u_d(p-m,2)+func(x(2),(n+1)*del_t+(time_counter-1),e)*del_t;
B(N-1,1)=u(n,N)+del_t*u_d(p-m,N)+func(x(N),(n+1)*del_t+(time_counter-1),e)*del_t;%
for j=2:N-2
    B(j,1)=u(n,j+1)+del_t*u_d((p-m),j+1)+func(x(j+1),(n+1)*del_t+(time_counter-1),e)*del_t; %
end
B;
u(n+1,2:N)=A\B;
end
u;
u_sol=u;
B;