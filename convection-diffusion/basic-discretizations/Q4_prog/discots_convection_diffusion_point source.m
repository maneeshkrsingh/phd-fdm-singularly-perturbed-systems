clc
clear all
format short

for ep_counter=1:5
    e=10^(-(ep_counter+1));
 for N_counter=1:5
     N=2^(N_counter+1);
     
      sig1=2.2*e*log(N);
      sig2=2.2*e*log(N);
      
A=zeros(N-1);
U=zeros(N+1,1);
D=zeros(N-1,1);
B=zeros(N-1,1);
x=zeros(1,N+1);
h=zeros(1,N);      
      
for i=1:(N/4)
    h(i)=4*(sig1)/N;
end
for i=(N/4)+1:N/2
    h(i)=4*(0.5-sig1)/N;
 end
for i=(N/2)+1:3*N/4
    h(i)=(4*sig2)/N;
end

for i=3*N/4+1:N
    h(i)=4*(0.5-sig2)/N;
end

 for i=1:N
    x(1)=0;
    x(i+1)=x(i)+h(i);
 end
 for l=1:2 
     if (l==2)
         N=2*N;
         sig1=2.2*e*log(N/2);
          sig2=2.2*e*log(N/2);
 
for i=1:(N/4)
    h(i)=4*(sig1)/N;
end
for i=(N/4)+1:N/2
    h(i)=4*(0.5-sig1)/N;
 end
for i=(N/2)+1:3*N/4
    h(i)=(4*sig2)/N;
end

for i=3*N/4+1:N
    h(i)=4*(0.5-sig2)/N;
end

 for i=1:N
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
     end

 
 
  U(1,1)=0;
  U(N+1,1)=0;

  %%%%%%%%%%%%%%%%%%%coefficints and sourse terms%%%%%%%%%%%%  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  for j=2:N
%     k=j-1;
%     f=(x(j+1)-x(j));
% %     h(j)
%     b=(x(j)-x(j-1));
% %     h(j-1)
%     c=(x(j+1)-x(j-1));
% 
%      if(x(j)<=0.5)
%         a1=1;
%      end
%     if(x(j)>.5)
%         a2=1;
%     end
%     %b1=-(x(j)*(1-x(j)));
%     b1=0;
%      if(x(j)<=0.5)
%         f1=2*(1+x(j)^2);
%      end
%     if(x(j)>.5)
%         f2=3*(1+x(j)^2);
%     end  
     
    
% a1=-1+x(j);
% a2=1+x(j);
 %   b1=0;
    if(k~=1 && k~=N-1)
        if(k<=N/2-1)
             A(k,k-1)=-e/(h(i)*h(i));
             A(k,k)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
             A(k,k+1)=-e/(h(i)*h(i+1))-1/h(i+1);
        end
        
        if(k>=N/2+1)
             A(k,k-1)=-e/(h(i)*h(i));
             A(k,k)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
             A(k,k+1)=-e/(h(i)*h(i+1))-1/h(i+1);
        end
       
        if(k==N/2)
             A(k,k-1)=-e/(h(i)*h(i));
             A(k,k)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
             A(k,k+1)=-e/(h(i)*h(i+1))-1/h(i+1);
        end
    end
    if(k==1)
              A(k,k)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
             A(k,k+1)=-e/(h(i)*h(i+1))-1/h(i+1);
    end
    if(k==N-1)
              A(k,k-1)=-e/(h(i)*h(i));
             A(k,k)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
    end        
        
   % end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  for i=2:N/2-1
%      A(i,i-1)=-2*e/((h(i)+h(i+1))*h(i))+a1(i)/h(i);
%      A(i,i)=2*e/(h(i)*h(i+1))-a1(i)/h(i)+b1(i);
%      A(i,i+1)=-2*e/((h(i)+h(i+1))*h(i+1));
%  end
%  
%  for i=N/2+1:N-2
%      A(i,i-1)=-2*e/((h(i)+h(i+1))*h(i));
%      A(i,i)=2*e/(h(i)*h(i+1))+a2(i)/h(i+1)+b2(i);
%      A(i,i+1)=-2*e/((h(i)+h(i+1))*h(i+1))-a2(i)/h(i);
%  end
%  
%   i=N/2;
%      A(i,i-1)=-1/h(i);
%      A(i,i)=(1/h(i)+1/h(i+1));
%      A(i,i+1)=-1/h(i+1);   
%      
%  i=1;
%      A(i,i)=2*e/(h(i)*h(i+1))-a1(i)/h(i)+b1(i);
%      A(i,i+1)=-2*e/((h(i)+h(i+1))*h(i+1));
%      
%  i=N-1;   
%      A(i,i-1)=-2*e/((h(i)+h(i+1))*h(i));
%      A(i,i)=2*e/(h(i)*h(i+1))+a2(i)/h(i+1)+b2(i);

%  for i=1:N-1
%      D(i)=f1(i+1);
%  end
 
%  for i=1:N/2
%      D(i)=f1(i+1);
% end 
%  for i=N/2+1:N-1
%      D(i)=f2(i+1);
%  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N-1
    j=k+1;
%         f=(x(j+1)-x(j));
%         b=(x(j)-x(j-1));
%         c=(x(j+1)-x(j-1));
     if(k~=1 && k~=N-1)
        if(k~=N/2)
             B(k,1)=x(k+1);
        end
        if(k==N/2)
            B(k,1)=x(k+1)+1/h(k+1);
        end
     end
     if(k==1)
         %a1=-(1+x(2)*(1-x(2)));
         B(k,1)=x(k+1);
     end
     if(k==N-1)
         %a2=(1+x(N)*(1-x(N)));
        B(k,1)=x(k+1);
     end
         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 U(2:N,1)=A\B;
if (l==1)
 U_pre=U;
else
    err=abs(U(1:2:N+1)-U_pre);
    max_err=max(err);
end

 end
  error(ep_counter,N_counter)=max_err;
 end

 end

for j=1:ep_counter
    for i=1:N_counter-1
         cnv(j,i)=log2(error(j,i)/error(j,i+1));
    end
end
plot(x,U)
error
cnv
