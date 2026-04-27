clc
%clear all
format short e
tic
%for ep_counter=1:3
    clear A
    clear B
    clear C
    clear D
%     e2=10^((ep_counter-1));
%     e1=10^((ep_counter));
    e=2^(-8);
   % e=2^(-2*(ep_counter+2));
    %e1=2^(-2*(ep_counter+10));
    %e2=2^(-2*(ep_counter+1));
    %e1=e2/2;
    %m1=exp(-1/e2);
    %m2=1-exp(-1/e2);
    %c1=1/(1-exp(-1/e1));
   % c2=1/(1-exp(-1/e2));

for N_counter =1:1
     %Nx1=16*2^(N_counter-1);
     %Nt1=10*4^(N_counter-1);
  
    Nx=512;
    Nt=256;
    
sig2=min(1/2,4.2*e*log(Nx));

%for l=1:1
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
q=zeros(Nx+1,1);
qp=zeros(Nx+1,1);
h=zeros(Nx,1);
t=zeros(Nt+1,1);
%c=zeros(Nx-1,1);
%%%%%%%%%%%%%%%%%%%%%%%


% construction of mesh points along  spatial variable and time variable
 
for i=1:(Nx/2)
      h(i)=2*(1-sig2)/Nx;
end

% for i=(Nx/2)+1:3*Nx/4
%     h(i)=(4*((sig2)-(sig1)))/Nx;
%  end
for i=(Nx/2)+1:Nx
    h(i)=(2*(sig2))/Nx;
end


 %%%%%%%%%%
for i=1:Nx-1
    x(1)=0;
    x(i+1)=x(i)+h(i);
    x(Nx+1)=1;
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
%     U1(i,1)=x(i);
%      U2(i,1)=x(i)*(1-x(i));
      U1(i,1)=0;
      U2(i,1)=0;
   % U2(i,1)=(c2*(1-exp(-x(i)/e2))-x(i));
end
 for j=1:Nt
    U1(1,j)=0;
    U1(Nx+1,j)=0;
    
    U2(1,j)=0;
    U2(Nx+1,j)=0;
 end
for j=1:Nt+1
    for i=1:Nx+1
       % f1(i,j)=(t(j)^3)*(x(i)^2)*(1-t(j))*(1-x(i));
        %f2(i,j)=(t(j)^3)*(x(i)^2)*(1-x(i))^2;
     
        f1(i,j)=(t(j)^3)*(1-t(j));
        f2(i,j)=(x(i)^2)*(1-x(i))^2;
        %f2(i,j)=((x(i)^2))*(t(j))*exp(-t(j));
    end
end

%%%%%%%%%%%%%%%%%%%% parameter %%%%%%%%%%%%%%
for i=1:Nx-1
if (7*h(i)>=2*e)
    q(i)=1/2;
else
    %q(i)=(h(i)-h(i+1))/3*h(i);
     q(i)=0;
end
    qp(i)=1-q(i);
end
   q(Nx)=0;
   qp(Nx)=1;

% for i=N/2+1:N
%     qm(i)=1/2;
% end

%%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 for j=2:Nt+1
for i=3:Nx-1
    A(i-1,i-2)=(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i-1)+q(i-2)*(-7*(2*h(i)+h(i-1))+h(i-1)*(h(i)+h(i-1))*(9+x(i-2))))+q(i-2);
    A(i-1,i-1)=q(i)*(9+x(i-2))+qp(i)*(9+x(i-1))-(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i)+q(i-1)*(-7*(2*h(i)+h(i-1))...
               +h(i)*(h(i)+h(i-1))*(9+x(i-1))))-(k/(h(i-1)*(h(i)+h(i-1))))*(-2*e+7*h(i-1)-14*h(i-1)*q(i-1));
    A(i-1,i)=(k/(h(i)*(h(i-1)+h(i))))*(-2*e+7*h(i)-14*h(i)*q(i));
end

i=2;
    A(i-1,i-1)=q(i)*(9+x(i-1))+qp(i)*(9+x(i))-(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i)+q(i-1)*(-7*(2*h(i)+h(i-1))...
               +h(i)*(h(i)+h(i-1))*(9+x(i-1))))-(k/(h(i-1)*(h(i)+h(i-1))))*(-2*e+7*h(i-1)-14*h(i-1)*q(i-1));
    A(i-1,i)=(k/(h(i)*(h(i-1)+h(i))))*(-2*e+7*h(i)-14*h(i)*q(i));
    
 i=Nx;   
    A(i-1,i-2)=(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i-1)+q(i-2)*(-7*(2*h(i)+h(i-1))+h(i-1)*(h(i)+h(i-1))*(9+x(i-2))))+q(i-2);
    A(i-1,i-1)=q(i)*(9+x(i-2))+qp(i)*(9+x(i-1))-(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i)+q(i-1)*(-7*(2*h(i)+h(i-1))...
               +h(i)*(h(i)+h(i-1))*(9+x(i-1))))-(k/(h(i-1)*(h(i)+h(i-1))))*(-2*e+7*h(i-1)-14*h(i-1)*q(i-1));
     
 for i=3:Nx-1
    D(i-1,i-2)=(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i-1)+q(i-2)*(-7*(2*h(i)+h(i-1))+h(i-1)*(h(i)+h(i-1))*(5+x(i-2))))+q(i-2);
    D(i-1,i-1)=q(i)*(9+x(i-2))+qp(i)*(5+x(i-1))-(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i)+q(i-1)*(-7*(2*h(i)+h(i-1))...
               +h(i)*(h(i)+h(i-1))*(5+x(i-1))))-(k/(h(i-1)*(h(i)+h(i-1))))*(-2*e+7*h(i-1)-14*h(i-1)*q(i-1));
    D(i-1,i)=(k/(h(i)*(h(i-1)+h(i))))*(-2*e+7*h(i)-14*h(i)*q(i));
end

i=2;
   D(i-1,i-1)=q(i)*(9+x(i-1))+qp(i)*(5+x(i-1))-(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i)+q(i-1)*(-7*(2*h(i)+h(i-1))...
               +h(i)*(h(i)+h(i-1))*(5+x(i-1))))-(k/(h(i-1)*(h(i)+h(i-1))))*(-2*e+7*h(i-1)-14*h(i-1)*q(i-1));
    D(i-1,i)=(k/(h(i)*(h(i-1)+h(i))))*(-2*e+7*h(i)-14*h(i)*q(i));
    
 i=Nx;   
    D(i-1,i-2)=(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i-1)+q(i-2)*(-7*(2*h(i)+h(i-1))+h(i-1)*(h(i)+h(i-1))*(5+x(i-2))))+q(i-2);
    D(i-1,i-1)=q(i)*(9+x(i-2))+qp(i)*(5+x(i-1))-(k/(h(i-1)*(h(i-1)+h(i))))*(-2*e-7*h(i)*qp(i)+q(i-1)*(-7*(2*h(i)+h(i-1))...
               +h(i)*(h(i)+h(i-1))*(5+x(i-1))))-(k/(h(i-1)*(h(i)+h(i-1))))*(-2*e+7*h(i-1)-14*h(i-1)*q(i-1));

  for i=3:Nx/2-1
     B(i-1,i-2)=-k*8*q(i-2);
     B(i-1,i-1)=-k*8*qp(i-1);
     B(i-1,i)=0;
  end
 
  
 i=2;
      B(i-1,i-1)=-k*8*qp(i-1);
     B(i-1,i)=0;
     
 i=Nx;   
     B(i-1,i-2)=-k*8*q(i-2);
     B(i-1,i-1)=-k*8*qp(i-1);
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% C matrix %%%%%%%%%  
     for i=3:Nx/2-1
     C(i-1,i-2)=-k*4*q(i-2);
     C(i-1,i-1)=-k*4*qp(i-1);
     C(i-1,i)=0;
     end
 
  
 i=2;
     C(i-1,i-1)=-k*4*qp(i-1);
     C(i-1,i)=0;
     
 i=Nx;   
     C(i-1,i-2)=-k*4*q(i-2);
     C(i-1,i-1)=-k*4*qp(i-1);
    
     M=[A,B;C,D];
  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %for j=2:Nt  
for i=3:Nx-1
             b1(i-1)=k*(q(i-1)*f1(i-1,j)+qp(i)*f1(i,j))+(q(i-1)*U1(i-1,j-1)+qp(i)*U1(i,j-1));
             b2(i-1)=k*(q(i-1)*f2(i-1,j)+qp(i)*f2(i,j))+(q(i-1)*U2(i-1,j-1)+qp(i)*U2(i,j-1));
end
     for i=(Nx/2)+2:Nx-1
         b1(i-1)=k*f1(i,j)+U1(i,j-1);
         b2(i-1)=k*f2(i,j)+U2(i,j-1);
     end
      b1(1)=k*(q(i-1)*f1(1,j)+qp(i)*f1(2,j))+(q(i-1)*U1(1,j-1)+qp(i)*U1(2,j-1));
      b2(1)=k*(q(i-1)*f2(1,j)+qp(i)*f2(2,j))+(q(i-1)*U2(1,j-1)+qp(i)*U2(2,j-1));
      b1(Nx-1)=qp(i)*f1(Nx+1,j);
      b2(Nx-1)=qp(i)*f2(Nx+1,j);
      b=[b1;b2];
         
    %  P=inv(M)*b;
      P=M\b;
      %%%%%%%%%%%%%%%%%%%%%%%% numerical solution

      for i=2:(Nx)
          U1(i,j)=P(i-1);
      end
   
      for i=Nx+1:2*(Nx)-1
          U2(i-Nx+1,j)=P(i-1);
      end
  end
%end
%      if(l==1)
%            F1=U1;
%           F2=U2;
%     end
%        if(l==2)
%           G1=U1;
%           G2=U2;
%        end
       
end



toc