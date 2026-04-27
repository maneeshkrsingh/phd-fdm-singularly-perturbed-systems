clc
clear all
format short e

for ep_counter=1:10
    clear A
    clear B
    clear C
    clear D
   e=2^(-2*(ep_counter));
     m1=exp(-1/e);
    m2=1-m1;
    c=1/(1-exp(-1/e));
    
for N_counter =1:4
    Nx=32*2^(N_counter-1);
    Nt=Nx/4;

% sig1=min(1/2,2*e*log(Nx));
% sig1;



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

% construction of mesh points along  spatial variable and time variable
  x(1)=0;
for i=1:(Nx/2)
    x(i+1)=-2*e*log(1-2*(1-(1/Nx))*(i/Nx));
end
% for i=(N/4)+1:(N/2)
%     x(i+1)=x((N/4)+1)-2*(e2-e1)*log(1-4*(1-(1/N))*(i/N-1/4));
% end
for i=(Nx/2)+1:Nx-1
   x(i+1)=1-((1-2*e*log(Nx))*2*(Nx-i)/Nx);
end
x(Nx+1)=1;

x;
for i=1:Nx
h(i)=x(i+1)-x(i);
end
% for i=1:(Nx/2)
%       h(i)=(2*sig1)/Nx;
% end
% 
% for i=(Nx/2)+1:Nx
%     h(i)=(2*(1-sig1))/Nx;
% end
% 
% 
%  %%%%%%%%%%
% for i=1:Nx
%     x(1)=0;
%     x(i+1)=x(i)+h(i);
% end
% x;
for j=1:Nt+1
    t(j)=(j-1)/Nt;
end
k=1/Nt;
t;
k;
%%%%%%%%%%%%%%%%%%%
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term and exact solution

 for i=1:Nx+1
    U1(i,1)=m1+m2*(1-x(i))-exp(-x(i)/e);
    U2(i,1)=c*(1-exp(-x(i)/e))-x(i)*exp(x(i)-1);
end
 for j=1:Nt
    U1(1,j)=0;
    U1(Nx+1,j)=0;
    
    U2(1,j)=0;
    U2(Nx+1,j)=0;
 end
for j=1:Nt+1
    for i=1:Nx+1
     
        f1(i,j)=exp(-t(j))*(exp(-x(i)/e)*(-(x(i)/e)-2*x(i)+c)+x(i)*exp(x(i)-1)+(2*m1*x(i)+m2*(3*x(i)-2*(x(i)^2)+1)-c));
        f2(i,j)=exp(-t(j))*(exp(-x(i)/e)*(-2*c*(x(i)/e)+x(i)-c*x(i))+exp(x(i)-1)*((1+2*x(i))*(1+x(i))+e*(x(i)+2)-(x(i)^2))+(x(i)*c-m1*x(i)-m2*x(i)*(1-x(i))));
end
end
%%%%%%%%%%%%%%%%%%
for j=1:Nt+1
    for i=1:Nx+1
       U1_exct(i,j)=exp(-t(j))*(m1+m2*(1-x(i))-exp(-x(i)/e));
       U2_exct(i,j)=exp(-t(j))*(c*(1-exp(-x(i)/e))-x(i)*exp(x(i)-1));
    end
end
%%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 % A
 for j=2:Nt+1
 for i=2:Nx-2
     A(i,i-1)=-2*e*k/((h(i+1)+h(i))*h(i));
     A(i,i)=(2*e*k/(h(i+1)*h(i)))+k*((1+x(i))/h(i+1))+k*(1+2*x(i))+1;
     A(i,i+1)=(-2*e*k/((h(i+1)+h(i))*h(i+1)))-k*(1+x(i))/h(i+1);
 end
 i=1;
      A(i,i)=(2*e*k/(h(i+1)*h(i)))+k*((1+x(i))/h(i+1))+k*(1+2*x(i))+1;
     A(i,i+1)=(-2*e*k/((h(i+1)+h(i))*h(i+1)))-k*(1+x(i))/h(i+1);
     
 i=Nx-1;   
    A(i,i-1)=-2*e*k/((h(i+1)+h(i))*h(i));
     A(i,i)=(2*e*k/(h(i+1)*h(i)))+k*((1+x(i))/h(i+1))+k*(1+2*x(i))+1;
%D
 for i=2:Nx-2
     D(i,i-1)=-2*e*k/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*e*k/(h(i+1)*h(i)))+k*((1+2*x(i))/h(i+1))+k*(1+(x(i)))+1;
     D(i,i+1)=(-2*e*k/((h(i+1)+h(i))*h(i+1)))-k*(1+2*x(i))/h(i+1);
 end
 i=1;
      D(i,i)=(2*e*k/(h(i+1)*h(i)))+k*((1+2*x(i))/h(i+1))+k*(1+(x(i)))+1;
     D(i,i+1)=(-2*e*k/((h(i+1)+h(i))*h(i+1)))-k*(1+2*x(i))/h(i+1);
     
 i=Nx-1;   
     D(i,i-1)=-2*e*k/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*e*k/(h(i+1)*h(i)))+k*((1+2*x(i))/h(i+1))+k*(1+(x(i)))+1;
     

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
     C(i,i)=-k*(x(i));
     C(i,i+1)=0;
 end
 i=1;
     C(i,i)=-k*(x(i));
     C(i,i+1)=0;
     
 i=Nx-1;   
     C(i,i-1)=0;
     C(i,i)=-k*(x(i));  
     

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
     
      P=M\b;
       %%%%%%%%%%%%%%%%%%%%%%%% numerical solution

      for i=2:(Nx)
          U1(i,j)=P(i-1);
      end
   
      for i=Nx+1:2*(Nx)-1
          U2(i-Nx+1,j)=P(i-1);
      end
 end

  %%%%%%%%%%%%%%%%%%%%%
   err1_normal=abs(U1-U1_exct);
    err2_normal=abs(U2-U2_exct);



err1_sup_normal=max(err1_normal);
max_err1_normal=max(err1_sup_normal);
error1_normal(ep_counter,N_counter)=max_err1_normal;

err2_sup_normal=max(err2_normal);
max_err2_normal=max(err2_sup_normal);
error2_normal(ep_counter,N_counter)=max_err2_normal;      
end
end

for j=1:ep_counter
    for i=1:N_counter-1
         conv1_normal(j,i)=log2(error1_normal(j,i)/error1_normal(j,i+1));
         conv2_normal(j,i)=log2(error2_normal(j,i)/error2_normal(j,i+1));
    end
end
error1_normal
error2_normal
conv1_normal
conv2_normal