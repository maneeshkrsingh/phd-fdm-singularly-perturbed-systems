clc
clear all
format long

%for ep_counter=1:10
    clear A
    clear B
    clear C
    clear D
    e2=2^(-5);
    e1=2^(-10);
    %e2=2^(-2*(ep_counter+2));
    %e1=2^(-2*(ep_counter+3));
     c1=1/(1-exp(-1/e1));
    c2=1/(1-exp(-1/e2));
 for N_counter=1:5
     N1=4^(N_counter+1);
     M1=N1/8;
     
      sig1=1.2*e1*log(N1);
      sig2=1.2*e2*log(N1);
for l=1:2
     N=N1*2^(l-1);
     M=M1*2^(l-1); 
      
A=zeros(N-1);
B=zeros(N-1);
C=zeros(N-1);
D=zeros(N-1);
U1=zeros(N+1,M+1);
U2=zeros(N+1,M+1);
d1=zeros(N-1,1);
d2=zeros(N-1,1);
M_mat=zeros(2*N-2,2*N-2);
x=zeros(N+1,1);
h=zeros(N,1);
t=zeros(M+1,1);
    %%%%%%%%%%%%%%%%%%%%%%%meshconstruction%%%%%%%%%%%%%%%%%%%  
for i=1:(N/4)
    h(i)=4*(0.5-sig2)/N;
end
for i=(N/4)+1:3*N/8
    h(i)=8*(sig2-sig1)/N;
 end
for i=(3*N/8)+1:N/2
    h(i)=(8*sig1)/N;
end
for i=(N/2)+1:5*N/8
    h(i)=(8*sig1)/N;
end
for i=5*N/8+1:6*N/8
    h(i)=8*(sig2-sig1)/N;
end
for i=(6*N/8)+1:N
    h(i)=4*(0.5-sig2)/N;
end
 for i=1:N
    x(1)=0;
    x(i+1)=x(i)+h(i);
 end
 for j=1:M+1
    t(j)=((j-1))/M;
end
del_t=1/M;

 
 %%%%%%%%%%%%%%%%%Boundary conditions and initai condition%%%%%%%%%%%%%%%%%%%%%%%%
   for i=1:N+1
     U1(i,1)=0;
     U2(i,1)=0;
    %U1(i,1)=(c1*(1-exp(-x(i)/e1))-x(i));
    % U2(i,1)=(c2*(1-exp(-x(i)/e2))-x(i));   
   end
 
  for j=1:M
    U1(1,j)=0;
    U1(N+1,j)=0;
    
    U2(1,j)=0;
    U2(N+1,j)=0;
 end
  %%%%%%%%%%%%%%%%%%%coefficints and sourse terms%%%%%%%%%%%%  
 
 
for i=2:M+1 
   for j=2:N
    k=j-1;
    f=(x(j+1)-x(j));%h(i+1)

    b=(x(j)-x(j-1));%h(i)

    c=(x(j+1)-x(j-1));%h(i+1)+h(i)

     if(x(j)<=0.5)
        a11=(1+x(j));
        a12=(1+x(j)^2);
     end
    if(x(j)>0.5)
        a21=-(2+x(j));
        a22=-(2+x(j)^2);
    end
   
    b11=4+x(j);
    b12=-2+x(j);
    b21=-4-x(j);
    b22=4+2*x(j);
  %%%%%%%%%%%%%%%%%%%%%A Matrix%%%%%%%%%%%%%  
    if(k~=1 && k~=N-1)
        if(k<=N/2-1)
             A(k,k-1)=del_t*(-2*e1/(c*b)-(a11/b));
             A(k,k)=del_t*(2*e1/(c*f)+2*e1/(c*b)+(a11/b)+b11)+1;
             A(k,k+1)=-del_t*2*e1/(f*c);
        end
        
        if(k>=N/2+1)
             A(k,k-1)=-del_t*2*e1/(c*b);
             A(k,k)=del_t*(2*e1/(c*f)+2*e1/(c*b)-(a21/f)+b11)+1;
             A(k,k+1)=del_t*(-2*e1/(c*f)+(a21/f));
        end
       
        if(k==N/2)
            A(k,k-1)=-1/b;
            A(k,k)=(1/b + 1/f);
            A(k,k+1)=-1/f;
        end
    end    
      if(k==1)
             A(k,k)=del_t*(2*e1/(c*f)+2*e1/(c*b)+(a11/b)+b11)+1;
             A(k,k+1)=-del_t*2*e1/(f*c);
     end
    if(k==N-1)
              A(k,k-1)=-del_t*2*e1/(c*b);
             A(k,k)=del_t*(2*e1/(c*f)+2*e1/(c*b)-(a21/f)+b11)+1;
    end           
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%D Matrix%%%%%%%%%%%%   
   if(k~=1 && k~=N-1)   
     if(k<=N/2-1)
             D(k,k-1)=-del_t*(2*e2/(c*b)+(a12/b));
             D(k,k)=del_t*(2*e2/(c*f)+2*e2/(c*b)+(a12/b)+b22)+1;
             D(k,k+1)=del_t*(-2*e2/(f*c));
     end
        
        if(k>=N/2+1)
             D(k,k-1)=-del_t*2*e2/(c*b);
             D(k,k)=del_t*(2*e2/(c*f)+2*e2/(c*b)-(a22/f)+b22)+1;
             D(k,k+1)=del_t*(-2*e2/(c*f)+(a22/f));
        end
       
        if(k==N/2)
            D(k,k-1)=-1/b;
            D(k,k)=(1/b + 1/f);
            D(k,k+1)=-1/f;
        end
    end
    if(k==1)
             D(k,k)=del_t*(2*e2/(c*f)+2*e2/(c*b)+(a12/b)+b22)+1;
             D(k,k+1)=del_t*(-2*e2/(f*c));
    end
    if(k==N-1)
             D(k,k-1)=-del_t*2*e2/(c*b);
             D(k,k)=del_t*(2*e2/(c*f)+2*e2/(c*b)-(a22/f)+b22)+1;
    end    
    %%%%%%%%%%%%%%%%%B matrix%%%%%%%%%%%%%%%%%%%%%
     if(k~=1 && k~=N-1)
             B(k,k-1)=0;
             B(k,k)=del_t*(b12);
             B(k,k+1)=0;
         
     end
    if(k==1)
             B(k,k)=del_t*(b12);
             B(k,k+1)=0;
    end
    if(k==N-1)
             B(k,k-1)=0;
             B(k,k)=del_t*(b12);
    end 
    %%%%%%%%%%%%%%%%%C matrix%%%%%%%%%%%%%%%%%%%%% 
     if(k~=1 && k~=N-1)
             C(k,k-1)=0;
             C(k,k)=del_t*(b21);
            C(k,k+1)=0;
         
     end
    if(k==1)
             C(k,k)=del_t*(b21);
             C(k,k+1)=0;
    end
    if(k==N-1)
             C(k,k-1)=0;
             C(k,k)=del_t*(b21);
    end 
   end   
 %%%%%%%%%%%%%%%%%end coeffiecient matrix%%%%%%%%%%%%%%%%%%%%%
% end
 M_mat=[A,B;C,D];

%%%%%%%%%%%%%%%%%%%Right hand side%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%first term%%%%%%%%%%%%%%
for k=1:N-1
    j=k+1;
        f=(x(j+1)-x(j));
        b=(x(j)-x(j-1));
        c=(x(j+1)-x(j-1));
     if(k~=1 && k~=N-1)
        if(k~=N/2)
             d1(k,1)=del_t*func1(x(k+1))+U1(k+1,i-1);
             d2(k,1)=del_t*func2(x(k+1))+U2(k+1,i-1);
        end
        if(k==N/2)
            d1(k,1)=0;
            d2(k,1)=0;
        end
     end
     if(k==1)
         %a11=(1+x(2));
         %a12=(1+x(2)^2);
         d1(k,1)=del_t*func1(x(k+1))+del_t*(2*e1/(c*b))*U1(1,i-1);
         d2(k,1)=del_t*func2(x(k+1))+del_t*(2*e2/(c*b))*U2(1,i-1);
     end
     if(k==N-1)
         %a21=-(2+x(N));
        % a22=-(2+x(N)^2);
        d1(k,1)=del_t*func1(x(k+1))+del_t*(2*e1/(c*f))*U1(N+1,i-1);
        d2(k,1)=del_t*func2(x(k+1))+del_t*(2*e2/(c*f))*U2(N+1,i-1);
     end
         
end

 d=[d1;d2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  P=M_mat\d;
       %%%%%%%%%%%%%%%%%%%%%%%% numerical solution

      for n=2:(N)
          U1(n,i)=P(n-1);
      end
   
      for n=N+1:2*(N)-1
          U2(n+1-N,i)=P(n-1);
      end
 % end
 end
if (l==1)
 F1=U1;
 F2=U2;
end
   if(l==2)
          G1=U1;
          G2=U2;
   end
       
end  
    
    err1=abs(G1(1:2:N+1,1:2:M+1)-F1);
     err_sup1=max(err1);
         max_err1=max(err_sup1);
         error1(N_counter)=max_err1;
      
      %error1(ep_counter,N_counter)=max_err1;
      
    err2=abs(G2(1:2:N+1,1:2:M+1)-F2);
    err_sup2=max(err2);
         max_err2=max(err_sup2);
         error2(N_counter)=max_err2;
      
       %error2(ep_counter,N_counter)=max_err2;
end

 %end
  

%for j=1:ep_counter
    for i=1:N_counter-1
         cnv1(i)=log2(error1(i)/error1(i+1));
         %cnv1(j,i)=log2(error1(j,i)/error1(j,i+1));
         cnv2(i)=log2(error2(i)/error2(i+1));
         %cnv2(j,i)=log2(error2(j,i)/error2(j,i+1));
    end
%end
% surf(x,t,U1')
 error1
 cnv1
 error2
 cnv2
% plot(x,U2(:,16),'r')
% hold on
% plot(x,U1(:,16))
% hold off

% figure
% subplot(1,2,1)
%  plot(x,U1(:,16))
%  %plot(x,y1)
% title('')
% 
% subplot(1,2,2)
%  plot(x,U2(:,16),'r')
% title('')
%toc
