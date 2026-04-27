
function [u_sol]= solution(e,N,x)   
A=zeros(N-1);
B=zeros(N-1,1);
u=zeros(1,N+1);
u_sol=zeros(1,N+1);

 %   u(1,:)=u_d(1,:);
%     
% for i=2:N+1
%     j=i-1;
%     h(j)=x(i)-x(i-1);
% end
% h;

    x;
    for j=2:N
    k=j-1;
    f=(x(j+1)-x(j));
%     h(j)
    b=(x(j)-x(j-1));
%     h(j-1)
    c=(x(j+1)-x(j-1));

     if(x(j)<=0.5)
        a1=-(1+x(j)*(1-x(j)));
     end
    if(x(j)>.5)
        a2=(1+x(j)*(1-x(j)));
    end
    %b1=-(x(j)*(1-x(j)));
    b1=1+x(j);
% a1=-1+x(j);
% a2=1+x(j);
 %   b1=0;
    if(k~=1 && k~=N-1)
        if(k<=N/2-1)
             A(k,k-1)=2*e/(c*b)-(a1/b);
             A(k,k)=-2*e/(c*f)-2*e/(c*b)+(a1/b)-b1;
             A(k,k+1)=2*e/(f*c);
        end
        
        if(k>=N/2+1)
             A(k,k-1)=2*e/(c*b);
             A(k,k)=-2*e/(c*f)-2*e/(c*b)-(a2/f)-b1;
             A(k,k+1)=2*e/(c*f)+(a2/f);
        end
       
        if(k==N/2)
            A(k,k-1)=1/b;
            A(k,k)=-(1/b + 1/f);
            A(k,k+1)=1/f;
        end
    end
    if(k==1)
             A(k,k)=-2*e/(c*f)-2*e/(c*b)+(a1/b)-b1;
             A(k,k+1)=2*e/(f*c);
    end
    if(k==N-1)
             A(k,k-1)=2*e/(c*b);
             A(k,k)=-2*e/(c*f)-2*e/(c*b)-(a2/f)-b1;
    end        
        
    end


u(1,1)=1; %Boundary condition
u(1,N+1)=0; %Boundary condition

%%%%%%%%%%%%%%
for k=1:N-1
    j=k+1;
        f=(x(j+1)-x(j));
        b=(x(j)-x(j-1));
        c=(x(j+1)-x(j-1));
     if(k~=1 && k~=N-1)
        if(k~=N/2)
             B(k,1)=func(x(k+1),e);
        end
        if(k==N/2)
            B(k,1)=0;
        end
     end
     if(k==1)
         a1=-(1+x(2)*(1-x(2)));
         B(k,1)=func(x(k+1),e)-(2*e/(c*b)-(a1/b))*u(1,1);
     end
     if(k==N-1)
         a2=(1+x(N)*(1-x(N)));
        B(k,1)=func(x(k+1),e)-(2*e/(c*f)+(a2/f))*u(1,N+1);
     end
         
end

%%%%%%%%%%%%%%%%%
% for j=1:N-1
%     if(j<=N/2)
%         B(j,1)=0.5*(u(n,j+1)-del_t*(func(x(j+1),n*del_t+(time_counter-1),e)))+0.5*(u(n,j)-del_t*(func(x(j),n*del_t+(time_counter-1),e)))-0.5*del_t*(u_d((p-m),j+1)+u_d((p-m),j));
%     else
%         B(j,1)=u(n,j+1)-del_t*func(x(j+1),n*del_t+(time_counter-1),e)-del_t*u_d((p-m),j+1);
%     end
% end
B;
u(1,2:N)=A\B;

u_sol=u;
