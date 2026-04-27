clc
clear all
format short

for ep_counter=1:5
    e=2^(-(ep_counter+3));
 for N_counter=1:2
     N1=16*2^(N_counter-1);
     
      sig1=4.2*e*log(N1);
      sig2=4.2*e*log(N1);
      
      
for l=1:2
      N=N1*2^(l-1);
      
A=zeros(N-1);
U=zeros(N+1,1);
%D=zeros(N-1,1);
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


%%%%%%%%%%%%%%%%%%%% Boundary conditions %%%%%%%%%%%%% 
 
  U(1,1)=0;
  U(N+1,1)=0;

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=2:N-2
             A(i,i-1)=-e/(h(i)*h(i));
             A(i,i)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
             A(i,i+1)=-e/(h(i)*h(i+1))-1/h(i+1);
 end
        
%   for i=N/2+2:N-2
%              A(i,i-1)=-e/(h(i)*h(i));
%              A(i,i)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
%              A(i,i+1)=-e/(h(i)*h(i+1))-1/h(i+1);
%   end
%        
%         if(i==N/2)
%              A(i,i-1)=-e/(h(i)*h(i));
%              A(i,i)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
%              A(i,k+1)=-e/(h(i)*h(i+1))-1/h(i+1);
%         end
   % end
    i=1;
              A(i,i)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
             A(i,i+1)=-e/(h(i)*h(i+1))-1/h(i+1);
    i=N-1;
              A(i,i-1)=-e/(h(i)*h(i));
             A(i,i)=(e/h(i))*(1/h(i+1)+1/h(i))+1/h(i+1);
   % end        
        
 %end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:N-1
   % j=k+1;

     if(k~=1 && k~=N-1)
        if(k~=N/2)
             B(k,1)=x(k+1);
        end
        if(k==N/2)
            B(k,1)=x(k+1)+1/h(k+1);
        end
     end
     if(k==1)
        
         B(k,1)=x(k+1);
     end
     if(k==N-1)
        
        B(k,1)=x(k+1);
     end
         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 U(2:N,1)=A\B;
 if(l==1)
           F=U;
         
 end
       if(l==2)
          G=U;
          
       end
 
end
    err=abs(F-G(1:2:N+1));
    err_sup=max(err);
         max_err=max(err_sup);
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
