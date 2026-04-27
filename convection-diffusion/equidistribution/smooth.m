function [v_sol]= smooth(e,N,p,x,del_t,v_d,time_counter)

A=zeros(N);
B=zeros(N,1);
v=zeros(p+1,N+1);
v_sol=zeros(p+1,N+1);
U=zeros(N+1,1);
D=zeros(N-1);

v_d;
v(1,:)=v_d(1,:);
v(:,1)=0; %boundary condition
for n=1:p
    for j=2:N+1
    k=j-1;
%     f=(x(j+1)-x(j));
%     h(j)
    b=(x(j)-x(j-1));
%     h(j-1)
%     c=(x(j+1)-x(j-1));
    if(k~=1)
        A(k,k)=1+(1+x(j)*(1-x(j)))*del_t/b;
        A(k,k-1)=-(1+x(j)*(1-x(j)))*del_t/b;
%         A(k,k+1)=-(2*e*del_t/(f*(b+f)));
    end
    if(k==1)
        A(k,k)=1+(1+x(j)*(1-x(j)))*del_t/b;
    end    
        
    end
  A;  
    m=n-1;
 n;
 p;
 del_t;
 v(n+1,1)=0; %Boundary condition
 
 B(1,1)=v(n,2)+del_t*v_d(p-m,2)+func(x(2),n*del_t+(time_counter-1),e)*del_t;%

 for j=2:N-1
     B(j,1)=v(n,j+1)+del_t*v_d(p-m,j+1)+func(x(j+1),n*del_t+(time_counter-1),e)*del_t;%
   
 end
B;
v(n+1,2:N+1)=inv(A)*B;

% for n=1:p
%     m=n-1;
%     for k=2:N+1
% %         if(n==1)
% %             v(n+1,k)=del_t*(-2*v_d(p+1-m,k))/exp(1)+v_d(n,k);
% %         else
%             v(n+1,k)=(del_t*v_d(p-m,k)+v(n,k)+(1+x(k)*(1-x(k)))*del_t*v(n+1,k-1)/(x(k)-x(k-1))+del_t*func(x(k),n*del_t+(time_counter-1),e))/(1+(1+x(k)*(1-x(k)))*del_t/(x(k)-x(k-1)));
% %         end
%     end
% end
end
v_sol=v;
v;
end
