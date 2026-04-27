clc
clear all
format short e

for ep_counter=1:2
    clear A
    clear B
    clear C
    clear D
    
    e=2^(-(9+ep_counter));
    
    m1=exp(-1/e);
    m2=1-m1;
    c=1/(1-exp(-1/e));
    

for N_counter =1:2
    Nx=64*2^(N_counter-1);
    %Nt=Nx;
    k=0.8/Nx;
    %k=1/Nt;
    Nt=1/k;
 
sig=min(1/2,4.2*e*log(Nx));

U1=zeros(Nx+1,Nt+1); 
U2=zeros(Nx+1,Nt+1); 
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
b=zeros(2*Nx-2,1);
P=zeros(2*Nx-2,1);
M=zeros(2*Nx-2,2*Nx-2);
x=zeros(Nx+1,1);
h=zeros(Nx,1);
t=zeros(Nt+1,1);

%%%%%%%%%%%%%%%%%%%%%%%
% construction of mesh points along  spatial variable and time variable
 
for i=1:(Nx/2)
      h(i)=(2*sig)/Nx;
end

for i=(Nx/2)+1:Nx
    h(i)=2*(1-sig)/Nx;
end


 %%%%%%%%%%
for i=1:Nx
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
x;
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
     
        f1(i,j)=exp(-t(j))*(exp(-x(i)/e)*(-(x(i)/e)-2*x(i)-c)+x(i)*exp(x(i)-1)+(2*m1*x(i)+m2*(3*x(i)-2*(x(i)^2)+1)-c));
        f2(i,j)=exp(-t(j))*(exp(-x(i)/e)*(-2*(x(i)/e)+x(i)-c*x(i))+exp(x(i)-1)*((1+2*x(i))*(e+x(i)+1)-x(i)^2)+(x(i)*c-m1*x(i)-m2*x(i)*(1-x(i))));
    end
end
%%%%%%%%%%%%%%%%%%
for j=1:Nt+1
    for i=1:Nx+1
       U1_exct(i,j)=exp(-t(j))*(m1+m2*(1-x(i))-exp(-x(i)/e));
       U2_exct(i,j)=exp(-t(j))*(c*(1-exp(-x(i)/e))-x(i)*exp(x(i)-1));
    end
end
%construction of finite diffrence matrix A B C D M
 % A
 for j=2:Nt+1
 for i=2:Nx/2-1
      A(i,i-1)=k*(-2*e/(h(i)*(h(i)+h(i+1)))+(1+x(i))/(h(i+1)+h(i)));
      A(i,i)=k*((2*e/(h(i+1)*h(i)))+(1+2*x(i)))+1;
      A(i,i+1)=k*((-2*e/(h(i+1)*(h(i)+h(i+1))))-(1+x(i))/(h(i+1)+h(i)));
 end
 i=1;
       A(i,i)=k*((2*e/(h(i+1)*h(i)))+(1+2*x(i)))+1;
      A(i,i+1)=k*((-2*e/(h(i+1)*(h(i)+h(i+1))))-(1+x(i))/(h(i+1)+h(i)));
 
 for i=(Nx/2):Nx-2
    A(i,i-1)=k*(-2*e/(h(i)*(h(i)+h(i+1))));
    A(i,i)=k*((2*e/(h(i+1)*h(i)))+((1+2*x(i))+(1+2*x(i+1)))/4+((1+x(i))+(1+x(i+1)))/(2*h(i+1)))+1/2;
    A(i,i+1)=k*(-2*e/(h(i+1)*(h(i)+h(i+1)))-(((1+x(i))+(1+x(i+1)))/(2*h(i+1)))+((1+2*x(i))+(1+2*x(i+1)))/4)+1/2;
end   
 i=Nx-1;      
     A(i,i-1)=k*(-2*e/(h(i)*(h(i)+h(i+1))));
    A(i,i)=k*((2*e/(h(i+1)*h(i)))+((1+2*x(i))+(1+2*x(i+1)))/4+((1+x(i))+(1+x(i+1)))/(2*h(i+1)))+1/2;
%D
for i=2:Nx/2-1
      D(i,i-1)=k*(-2*e/(h(i)*(h(i)+h(i+1)))+(1+2*x(i))/(h(i+1)+h(i)));
      D(i,i)=k*((2*e/(h(i+1)*h(i)))+(1+x(i)))+1;
      D(i,i+1)=k*((-2*e/(h(i+1)*(h(i)+h(i+1))))-(1+2*x(i))/(h(i+1)+h(i)));
end
 i=1;
      D(i,i)=k*((2*e/(h(i+1)*h(i)))+(1+x(i)))+1;
      D(i,i+1)=k*((-2*e/(h(i+1)*(h(i)+h(i+1))))-(1+2*x(i))/(h(i+1)+h(i)));
 
 for i=(Nx/2):Nx-2
    D(i,i-1)=k*(-2*e/(h(i)*(h(i)+h(i+1))));
    D(i,i)=k*((2*e/(h(i+1)*h(i)))+((1+x(i))+(1+x(i+1)))/4+((1+2*x(i))+(1+2*x(i+1)))/(2*h(i+1)))+1/2;
    D(i,i+1)=k*(-2*e/(h(i+1)*(h(i)+h(i+1)))-((1+2*x(i))+(1+2*x(i+1))/2*h(i+1))+((1+x(i))+(1+x(i+1)))/4)+1/2;
end   
 i=Nx-1;      
   D(i,i-1)=k*(-2*e/(h(i)*(h(i)+h(i+1))));
    D(i,i)=k*((2*e/(h(i+1)*h(i)))+((1+x(i))+(1+x(i+1)))/4+((1+2*x(i))+(1+2*x(i+1)))/(2*h(i+1)))+1/2;

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
     C(i,i)=-(0.5)*k*(x(i)+x(i+1));
     C(i,i+1)=0;
   end
 i=1;
     C(i,i)=-(0.5)*k*(x(i)+x(i+1));
     C(i,i+1)=0;
     
 i=Nx-1;   
     C(i,i-1)=0;
     C(i,i)=-(0.5)*k*(x(i)+x(i+1));  
     

   M=[A,B;C,D];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %    for j=2:Nt  
     for i=1:(Nx/2)-1
             b1(i)=(k*f1(i+1,j)+U1(i+1,j-1));
              b2(i)=(k*f2(i+1,j)+U2(i+1,j-1));
     end
     for i=(Nx/2):Nx-1
        b1(i)=0.5*k*(f1(i+1,j)+f1(i,j))+0.5*(U1(i+1,j-1)+U1(i,j-1));
        b2(i)=0.5*k*(f2(i+1,j)+f2(i,j))+0.5*(U2(i+1,j-1)+U2(i,j-1));
    end
%       b1(1)=(k*f1(1,j)+U1(2,j-1));
%       b2(1)=(k*f2(1,j)+U2(2,j-1));
%        b1(Nx-1)=0.5*k*(f1(i+1,j)+f1(i,j))+0.5*(U1(i+1,j-1)+U1(i,j-1));
%         b2(Nx-1)=0.5*k*(f2(i+1,j)+f2(i,j))+0.5*(U2(i+1,j-1)+U2(i,j-1));
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
  
%   %%%%%%%%%%%%%%%%%%%%%
   err1=abs(U1-U1_exct);
   err1_sup=max(err1);
   max_err1=max(err1_sup);
   error1(ep_counter,N_counter)=max_err1;

 err2=abs(U2-U2_exct);
 err2_sup=max(err2);
 max_err2=max(err2_sup);
 error2(ep_counter,N_counter)=max_err2;      





end
end

for j=1:ep_counter
    for i=1:N_counter-1
         conv1(j,i)=log2(error1(j,i)/error1(j,i+1));
         conv2(j,i)=log2(error2(j,i)/error2(j,i+1));
    end
end
error1
error2
conv1
conv2