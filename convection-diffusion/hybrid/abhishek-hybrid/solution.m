
function [u_sol]= solution(e,N,p,x,del_t,u_d,time_counter)   
A=zeros(N-1,N-1);
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
        if(j<=N/2+1)
             A(k,k-1)=del_t*(-2*e/(b*c)-(((1+x(j-1)*(1-x(j-1)))+(1+x(j)*(1-x(j))))/(2*b)))+1/2;
             A(k,k)=del_t*(2*e/(b*f)+(((1+x(j-1)*(1-x(j-1)))+(1+x(j)*(1-x(j))))/(2*b)))+1/2;
             A(k,k+1)=del_t*(-2*e/(f*c));
        else
             A(k,k-1)=del_t*(-2*e/(b*c)-(1+x(j)*(1-x(j)))/c);
             A(k,k)=del_t*(2*e/(b*f))+1;
             A(k,k+1)=del_t*(-2*e/(f*c)+(1+x(j)*(1-x(j)))/c);
        end
    end
    if(k==1)
        A(k,k)=del_t*(2*e/(b*f)+(((1+x(j-1)*(1-x(j-1)))+(1+x(j)*(1-x(j))))/(2*b)))+1/2;
        A(k,k+1)=del_t*(-2*e/(f*c));
    end
    if(k==N-1)
        A(k,k)=del_t*(2*e/(b*f))+1;
        A(k,k-1)=del_t*(-2*e/(b*c)-(1+x(j)*(1-x(j)))/c);
    end        
        
    end
m=n-1;
A;
u(n+1,1)=0; %Boundary condition
u(n+1,N+1)=0; %Boundary condition

for j=1:N-1
    if(j<=N/2)
        B(j,1)=0.5*(u(n,j+1)+del_t*(func(x(j+1),n*del_t+(time_counter-1),e)))+0.5*(u(n,j)+del_t*(func(x(j),n*del_t+(time_counter-1),e)))+0.5*del_t*(u_d((p-m),j+1)+u_d((p-m),j));
    else
        B(j,1)=u(n,j+1)+del_t*func(x(j+1),n*del_t+(time_counter-1),e)+del_t*u_d((p-m),j+1);
    end
end
B;
u(n+1,2:N)=A\B;
end
u;
u_sol=u;
B;