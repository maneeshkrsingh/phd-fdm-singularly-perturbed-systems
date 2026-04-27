
function [u_sol]= solution(e,N,x)   
A=zeros(N-1);
B=zeros(N-1,1);
u=zeros(1,N+1);
u_sol=zeros(1,N+1);


    x;
    for j=2:N
    k=j-1;
    f=(x(j+1)-x(j));
    b=(x(j)-x(j-1));
    c=(x(j+1)-x(j-1));
    a1=1+x(j);
     b1=x(j);
    if(k~=1 && k~=N-1)
        if(k<=N/2-1)
             A(k,k-1)=2*e/(c*b);
             A(k,k)=(-2*e/(c*f))-(2*e/(c*b))-(a1/f)-b1;
             A(k,k+1)=(2*e/(f*c))+(a1/f);
        end
        
        if(k>=N/2+1)
             A(k,k-1)=2*e/(c*b);
             A(k,k)=(-2*e/(c*f))-(2*e/(c*b))-(a1/f)-b1;
             A(k,k+1)=(2*e/(c*f))+(a1/f);
        end
       
        if(k==N/2)
            A(k,k-1)=1/b;
            A(k,k)=-((1/b) + (1/f));
            A(k,k+1)=1/f;
        end
    end
    if(k==1)
             A(k,k)=-(2*e/(c*f))-(2*e/(c*b))-(a1/f)-b1;
             A(k,k+1)=(2*e/(f*c))+(a1/f);
    end
    if(k==N-1)
             A(k,k-1)=(2*e/(c*b));
             A(k,k)=-(2*e/(c*f))-(2*e/(c*b))-(a1/f)-b1;
    end        
        
    end


u(1,1)=-1; %Boundary condition
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
 %B(k,1)=-9;
        end
        if(k==N/2)
            B(k,1)=0;
        end
     end
     if(k==1)
         B(k,1)=func(x(k+1),e)-(2*e/(c*b))*u(1,1);
     end
     if(k==N-1)
        B(k,1)=func(x(k+1),e)-((2*e/(c*f))+(a1/f))*u(1,N+1);
     end
         
end

u(1,2:N)=A\B;

u_sol=u;
