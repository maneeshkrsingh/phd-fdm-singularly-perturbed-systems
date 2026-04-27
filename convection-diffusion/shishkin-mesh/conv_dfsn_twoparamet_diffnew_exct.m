clc
clear all
format short e

for ep_counter=1:10
    clear A
    clear B
    clear C
    clear D
    e2=2^(-2*(ep_counter-1));
    e1=2^(-2*(ep_counter+1));
%     m1=exp(-1/e2);
%     m2=1-m1;
    c2=1/(1-exp(-1/e2));
    c1=1/(1-exp(-1/e1));

for N_counter =1:3
    Nx=32*2^(N_counter-1);
    Nt=Nx/2;
 
sig2=min(1/2,4*e2*log(Nx));
sig1=min(1/4,4*e1*log(Nx));
sig2;
sig1;
% for l=1:2
%      Nx=Nx1*2^(l-1);
%     Nt=Nt1*2^(l-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U=zeors(2*Nx+2,Nt+1);
% U_num=zeros(2*Nx-2,Nt+1);
U1=zeros(Nx+1,Nt+1); %U1(2:Nx,:)=U_num(1:Nx-1,:);
%U_1=zeros(Nx-1,Nt);
U2=zeros(Nx+1,Nt+1); %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:); 
%U_2=zeros(Nx-1,Nt);
U1_exct=zeros(Nx+1,Nt+1);
U2_exct=zeros(Nx+1,Nt+1);
f1=zeros(Nx-1,1);
f2=zeros(Nx-1,1);
A=zeros(Nx-1);
B=zeros(Nx-1);
C=zeros(Nx-1);
D=zeros(Nx-1);
b1=zeros(Nx-1,1);
b2=zeros(Nx-1,1);
d1=zeros(Nx-1,1);
d2=zeros(Nx-1,1);
b=zeros(2*Nx-2,1);
P=zeros(2*Nx-2,1);
M=zeros(2*Nx-2,2*Nx-2);
x=zeros(Nx+1,1);
h=zeros(Nx,1);
t=zeros(Nt+1,1);
%c=zeros(Nx-1,1);
%%%%%%%%%%%%%%%%%%%%%%%
% construction of mesh points along  spatial variable and time variable
 
for i=1:(Nx/4)
      h(i)=(4*sig1)/Nx;
end

for i=(Nx/4)+1:Nx/2
    h(i)=(4*((sig2)-(sig1)))/Nx;
 end
for i=(Nx/2)+1:Nx
    h(i)=(2*(1-sig2))/Nx;
end


 %%%%%%%%%%
for i=1:Nx
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
x;
for j=1:Nt+1
    t(j)=((j-1))/Nt;
end
k=1/Nt;
t;
k;
%%%%%%%%%%%%%%%%%%%
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term and exact solution

 for i=1:Nx+1
   % U1(i,1)=(m1+m2*(1-x(i))-exp(-(x(i))/e2));
    U1(i,1)=(c1*(1-exp(-(x(i))/e1))+c2*(1-exp(-(x(i))/e2))-2*sin((pi/2)*x(i)));
    U2(i,1)=(c2*(1-exp(-(x(i))/e2))-x(i)*exp(x(i)-1));
end
 for j=1:Nt
    U1(1,j)=0;
    U1(Nx+1,j)=0;
    
    U2(1,j)=0;
    U2(Nx+1,j)=0;
 end
for j=1:Nt+1
    for i=1:Nx+1
     
        f1(i,j)=exp(-t(j))*(exp(-(x(i))/e1)*(-c1)+exp(-(x(i))/e2)*(-(c2/e2)+(c1*e1/(e2)^2))+(x(i)*exp(x(i)-1))+(c1+2-2*sin((pi/2)*x(i))));
        
         f2(i,j)=exp(-t(j))*(exp(-(x(i))/e1)*(c1)+exp(-(x(i))/e2)*(-2*c2-(c2/e2))+4*exp(x(i)-1)+(-c1+2*c2+2*sin((pi/2)*x(i))));
end
end
%%%%%%%%%%%%%%%%%%
for j=1:Nt+1
    for i=1:Nx+1
       U1_exct(i,j)=exp(-t(j))*(c1*(1-exp(-(x(i))/e1))+c2*(1-exp(-(x(i))/e2))-2*sin((pi/2)*x(i)));
       U2_exct(i,j)=exp(-t(j))*(c2*(1-exp(-(x(i))/e2))-x(i)*exp(x(i)-1));
       %U2_exct(i,j)=exp(-t(j))*(c1+c2*x(i)-exp(-(1-x(i))/e2));
    end
end
%%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 % A
 for j=2:Nt+1
 for i=2:Nx-2
     A(i,i-1)=-2*e1*k/((h(i+1)+h(i))*h(i));
     A(i,i)=(2*e1*k/(h(i+1)*h(i)))+k/(h(i+1))+2*k+1;
     A(i,i+1)=(-2*e1*k/((h(i+1)+h(i))*h(i+1)))-k/h(i+1);
 end
 i=1;
      A(i,i)=(2*e1*k/(h(i+1)*h(i)))+k/(h(i+1))+2*k+1;
     A(i,i+1)=(-2*e1*k/((h(i+1)+h(i))*h(i+1)))-k/h(i+1);
     
 i=Nx-1;   
   A(i,i-1)=-2*e1*k/((h(i+1)+h(i))*h(i));
    A(i,i)=(2*e1*k/(h(i+1)*h(i)))+k/(h(i+1))+2*k+1;
%D
 for i=2:Nx-2
     D(i,i-1)=-2*e2*k/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*e2*k/(h(i+1)*h(i)))+2*k/h(i+1)+4*k+1;
     D(i,i+1)=(-2*e2*k/((h(i+1)+h(i))*h(i+1)))-2*k/h(i+1);
 end
 i=1;
      D(i,i)=(2*e2*k/(h(i+1)*h(i)))+2*k/h(i+1)+4*k+1;
     D(i,i+1)=(-2*e2*k/((h(i+1)+h(i))*h(i+1)))-2*k/h(i+1);
     
 i=Nx-1;   
   D(i,i-1)=-2*e2*k/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*e2*k/(h(i+1)*h(i)))+2*k/h(i+1)+4*k+1;
     

  for i=2:Nx-2
     B(i,i-1)=0;
     B(i,i)=-k;
     B(i,i+1)=0;
 end
 i=1;
     B(i,i)=-k;
     B(i,i+1)=0;
     
 i=Nx-1;   
     B(i,i-1)=0;
     B(i,i)=-k;
     
   for i=2:Nx-2
     C(i,i-1)=0;
     C(i,i)=-k;
     C(i,i+1)=0;
 end
 i=1;
     C(i,i)=-k;
     C(i,i+1)=0;
     
 i=Nx-1;   
     C(i,i-1)=0;
     C(i,i)=-k;  
     

   M=[A,B;C,D];
  
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %    for j=2:Nt  
     for i=2:Nx-2
             b1(i)=(k*f1(i,j)+U1(i+1,j-1));
              b2(i)=(k*f2(i,j)+U2(i+1,j-1));
     end
%      
      b1(1)=(k*f1(1,j)+U1(2,j-1));
      b2(1)=(k*f2(1,j)+U2(2,j-1));
      b1(Nx-1)=k*f1(Nx-1,j)+U1(Nx,j-1);
      b2(Nx-1)=k*f2(Nx-1,j)+U2(Nx,j-1);
      b=[b1;b2];
     
      P=inv(M)*b;
       %%%%%%%%%%%%%%%%%%%%%%%% numerical solution

      for i=2:(Nx)
          U1(i,j)=P(i-1);
      end
   
      for i=Nx+1:2*(Nx)-1
          U2(i-Nx+1,j)=P(i-1);
      end
 end
  
%      if(l==1)
%            F1=U1;
%           F2=U2;
%     end
%        if(l==2)
%           G1=U1;
%           G2=U2;
%        end
%        
%   end
  %%%%%%%%%%%%%%%%%%%%%
 
   err1=abs(U1-U1_exct);
  %err1=abs(F1-G1(1:2:Nx+1,1:2:Nt+1));
    err_sup1=max(err1);
         max_err1=max(err_sup1);
       error1(ep_counter,N_counter)=max_err1;
          
           err2=abs(U2-U2_exct);
        %err2=abs(F2-G2(1:2:Nx+1,1:2:Nt+1));
        err_sup2=max(err2);
         max_err2=max(err_sup2);
       error2(ep_counter,N_counter)=max_err2;
%           
end
end
%end
for j=1:ep_counter
    for i=1:N_counter-1
         cnv_rt1(j,i)=log2(error1(j,i)/error1(j,i+1));
    end
end
for j=1:ep_counter
    for i=1:N_counter-1
         cnv_rt2(j,i)=log2(error2(j,i)/error2(j,i+1));
    end
end

cnv_rt1
cnv_rt2
error1
error2