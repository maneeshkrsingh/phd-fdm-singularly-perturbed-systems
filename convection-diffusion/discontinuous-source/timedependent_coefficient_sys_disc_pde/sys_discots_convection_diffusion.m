clc
clear long
format short

for ep_counter=1:10
    clear A
    clear B
    clear C
    clear D
    e2=2^(-2*(ep_counter+1));
    e1=2^(-2*(ep_counter+1));
 for N_counter=1:4
     N=4^(N_counter+1);
    
     
      sig1=2.2*e1*log(N);
      sig2=2.2*e2*log(N);
      
A=zeros(N-1);
B=zeros(N-1);
C=zeros(N-1);
D=zeros(N-1);
U1=zeros(N+1,1);
U2=zeros(N+1,1);
d1=zeros(N-1,1);
d2=zeros(N-1,1);
M=zeros(2*N-2,2*N-2);
x=zeros(1,N+1);
h=zeros(1,N);

    %%%%%%%%%%%%%%%%%%%%%%%meshconstruction%%%%%%%%%%%%%%%%%%%  
% for i=1:(N/4)
%     h(i)=4*(0.5-sig2)/N;
% end
% for i=(N/4)+1:3*N/8
%     h(i)=8*(sig2-sig1)/N;
%  end
% for i=(3*N/8)+1:N/2
%     h(i)=(8*sig1)/N;
% end
% for i=(N/2)+1:5*N/8
%     h(i)=(8*sig1)/N;
% end
% for i=5*N/8+1:6*N/8
%     h(i)=8*(sig2-sig1)/N;
% end
% for i=(6*N/8)+1:N
%     h(i)=4*(0.5-sig2)/N;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(N/4)
    h(i)=4*(0.5-sig1)/N;
end
for i=(N/4)+1:N/2
    h(i)=(4*sig1)/N;
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

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%doubling the mesh%%%%%%%%%%
 for l=1:2 
     if (l==2)
         N=2*N;
          sig1=2.2*e1*log(N/2);
          sig2=2.2*e2*log(N/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% for i=1:(N/4)
%     h(i)=4*(0.5-sig2)/N;
% end
% for i=(N/4)+1:3*N/8
%     h(i)=8*(sig2-sig1)/N;
%  end
% for i=(3*N/8)+1:N/2
%     h(i)=(8*sig1)/N;
% end
% for i=(N/2)+1:5*N/8
%     h(i)=(8*sig1)/N;
% end
% for i=5*N/8+1:6*N/8
%     h(i)=8*(sig2-sig1)/N;
% end
% for i=(6*N/8)+1:N
%     h(i)=4*(0.5-sig2)/N;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(N/4)
    h(i)=4*(0.5-sig1)/N;
end
for i=(N/4)+1:N/2
    h(i)=(4*sig1)/N;
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

 
 %%%%%%%%%%%%%%%%%Boundary conditions%%%%%%%%%%%%%%%%%%%%%%%%
  U1(1,1)=1;
  U1(N+1,1)=1;
  U2(1,1)=1;
  U2(N+1,1)=1;

  %%%%%%%%%%%%%%%%%%%coefficints and sourse terms%%%%%%%%%%%%  

 for j=2:N
    k=j-1;
    f=(x(j+1)-x(j));% h(j+1)
%     h(j)
    b=(x(j)-x(j-1));% h(j)
%     h(j-1)
    c=(x(j+1)-x(j-1));% h(j+1)+h(j)

     if(x(j)<=0.5)
        a11=-(1+x(j));
        a12=-(1+x(j)^2);
     end
    if(x(j)>.5)
        a21=(2+x(j));
        a22=(2+x(j)^2);
    end
   
    b11=4+x(j);
    b12=-2-x(j);
    b21=-4-x(j);
    b22=4+2*x(j);

 
 %%%%%%%%%%%%%%%%%%%%%A Matrix%%%%%%%%%%%%%
    if(k~=1 && k~=N-1)
        if(k<=N/2-1)
             A(k,k-1)=2*e1/(c*b)-(a11/b);
             A(k,k)=-2*e1/(c*f)-2*e1/(c*b)+(a11/b)-b11;
             A(k,k+1)=2*e1/(f*c);
        end
        
        if(k>=N/2+1)
             A(k,k-1)=2*e1/(c*b);
             A(k,k)=-2*e1/(c*f)-2*e1/(c*b)-(a21/f)-b11;
             A(k,k+1)=2*e1/(c*f)+(a21/f);
        end
       
        if(k==N/2)
            A(k,k-1)=1/b;
            A(k,k)=-(1/b + 1/f);
            A(k,k+1)=1/f;
        end
    end    
      if(k==1)
             A(k,k)=-2*e1/(c*f)-2*e1/(c*b)+(a11/b)-b11;
             A(k,k+1)=2*e1/(f*c);
     end
    if(k==N-1)
             A(k,k-1)=2*e1/(c*b);
             A(k,k)=-2*e1/(c*f)-2*e1/(c*b)-(a21/f)-b11;
    end           
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%D Matrix%%%%%%%%%%%%   
   if(k~=1 && k~=N-1)   
     if(k<=N/2-1)
             D(k,k-1)=2*e2/(c*b)-(a12/b);
             D(k,k)=-2*e2/(c*f)-2*e2/(c*b)+(a12/b)-b22;
             D(k,k+1)=2*e2/(f*c);
     end
        
        if(k>=N/2+1)
             D(k,k-1)=2*e2/(c*b);
             D(k,k)=-2*e2/(c*f)-2*e2/(c*b)-(a22/f)-b22;
             D(k,k+1)=2*e2/(c*f)+(a22/f);
        end
       
        if(k==N/2)
            D(k,k-1)=1/b;
            D(k,k)=-(1/b + 1/f);
            D(k,k+1)=1/f;
        end
    end
    if(k==1)
             D(k,k)=-2*e2/(c*f)-2*e2/(c*b)+(a12/b)-b22;
             D(k,k+1)=2*e2/(f*c);
    end
    if(k==N-1)
             D(k,k-1)=2*e2/(c*b);
             D(k,k)=-2*e2/(c*f)-2*e2/(c*b)-(a22/f)-b22;
    end    
    %%%%%%%%%%%%%%%%%B matrix%%%%%%%%%%%%%%%%%%%%%
     if(k~=1 && k~=N-1)
             B(k,k-1)=0;
             B(k,k)=-b12;
             B(k,k+1)=0;
         
     end
    if(k==1)
             B(k,k)=-b12;
             B(k,k+1)=0;
    end
    if(k==N-1)
             B(k,k-1)=0;
             B(k,k)=-b12;
    end 
    %%%%%%%%%%%%%%%%%C matrix%%%%%%%%%%%%%%%%%%%%% 
     if(k~=1 && k~=N-1)
             C(k,k-1)=0;
             C(k,k)=-b21;
            C(k,k+1)=0;
         
     end
    if(k==1)
             C(k,k)=-b21;
             C(k,k+1)=0;
    end
    if(k==N-1)
             C(k,k-1)=0;
             C(k,k)=-b21;
    end 
   
 %%%%%%%%%%%%%%%%%end coeffiecient matrix%%%%%%%%%%%%%%%%%%%%%
  end
 M=[A,B;C,D];

%%%%%%%%%%%%%%%%%%%Right hand side%%%%%%%%%%
for k=1:N-1
    j=k+1;
        f=(x(j+1)-x(j));
        b=(x(j)-x(j-1));
        c=(x(j+1)-x(j-1));
     if(k~=1 && k~=N-1)
        if(k~=N/2)
             d1(k,1)=func1(x(k+1));
             d2(k,1)=func2(x(k+1));
        end
        if(k==N/2)
            d1(k,1)=0;
            d2(k,1)=0;
        end
     end
     if(k==1)
         a11=-(1+x(2));
         a12=-(1+x(2)^2);
         d1(k,1)=func1(x(k+1))-(2*e1/(c*b)-(a11/b))*U1(1,1);
         d2(k,1)=func2(x(k+1))-(2*e2/(c*b)-(a12/b))*U2(1,1);
     end
     if(k==N-1)
         a21=(2+x(N));
         a22=(2+x(N)^2);
        d1(k,1)=func1(x(k+1))-(2*e1/(c*f)+(a21/f))*U1(N+1,1);
        d2(k,1)=func2(x(k+1))-(2*e2/(c*f)+(a22/f))*U2(N+1,1);
     end
         
end

 d=[d1;d2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  P=M\d;
       %%%%%%%%%%%%%%%%%%%%%%%% numerical solution

      for i=2:(N)
          U1(i,1)=P(i-1);
      end
   
      for i=N+1:2*(N)-1
          U2(i+1-N,1)=P(i-1);
      end


if (l==1)
 U1_pre=U1;
 U2_pre=U2;
else
    err1=abs(U1(1:2:N+1)-U1_pre);
    err2=abs(U2(1:2:N+1)-U2_pre);
    max_err1=max(err1);
    max_err2=max(err2);
end

 end
  error1(ep_counter,N_counter)=max_err1;
  error2(ep_counter,N_counter)=max_err2;
 end

 end

for j=1:ep_counter
    for i=1:N_counter-1
         cnv1(j,i)=log2(error1(j,i)/error1(j,i+1));
         cnv2(j,i)=log2(error2(j,i)/error2(j,i+1));
    end
end
plot(x,U1)
error1
cnv1
error2
cnv2
