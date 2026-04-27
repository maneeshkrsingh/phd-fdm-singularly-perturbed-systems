clc
clear all
format short e

%for ep_counter=1:10
    clear A
    clear B
    clear C
    clear D
     e2=2^(-20);
    e1=2^(-26);
    %e2=2^(-2*(ep_counter+3));
    %e1=2^(-2*(ep_counter+8));
    %e2=2^(-2*(ep_counter+1));
    %e1=e2/2;
    %m1=exp(-1/e2);
    %m2=1-exp(-1/e2);
    c1=1/(1-exp(-1/e1));
    c2=1/(1-exp(-1/e2));

for N_counter =1:2
    Nx1=16*2^(N_counter-1);
   % Nt1=10*Nx1/8;
    Nt1=(Nx1)^2;
 
sig2=min(1/2,1.4*e2*log(Nx1));
sig1=min(sig2/2,1.2*(e1)*log(Nx1));

% sig2=0.0024;
% sig1=0.0014;
for l=1:2
     Nx=Nx1*2^(l-1);
    Nt=Nt1*2^(l-1);

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
 
for i=1:(Nx/2)
      h(i)=2*(1-sig2)/Nx;
end

for i=(Nx/2)+1:3*Nx/4
    h(i)=(4*((sig2)-(sig1)))/Nx;
 end
for i=(3*Nx/4)+1:Nx
    h(i)=(4*(sig1))/Nx;
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
     U1(i,1)=(c1*(1-exp(-(1-x(i))/e1))-(1-x(i)));
      U2(i,1)=(c2*(1-exp(-(1-x(i))/e2))-(1-x(i)));
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
     
        f1(i,j)=((x(i)^2))*(t(j))*exp(-t(j));
        f2(i,j)=((x(i)^3)-x(i)+1)*(t(j)^2)*exp(-4*t(j));
        %f2(i,j)=((x(i)^2))*(t(j))*exp(-t(j));
         %f1(i,j)=1;
        %f2(i,j)=1;
    end
end

%%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 for j=2:Nt+1
for i=3:Nx/2+1
    A(i-1,i-2)=k*(-(2*e1/(h(i-1)*(h(i-1)+h(i))))-(((1+x(i))+(1+x(i-1)))/(2*h(i-1)))+((2+x(i))+(2+x(i-1)))/2)+(1/2);
    A(i-1,i-1)=k*((2*e1/(h(i-1)*h(i)))+(((1+x(i))+(1+x(i-1)))/(2*h(i-1)))+((2+x(i))+(2+x(i-1)))/2)+(1/2);
    A(i-1,i)=-k*(2*e1/(h(i)*(h(i-1)+h(i))));
end
for i=(Nx/2)+2:Nx-1
    A(i-1,i-2)=k*(-(2*e1/(h(i-1)*(h(i-1)+h(i))))-((1+x(i))/(h(i-1)+h(i))));
    A(i-1,i-1)=k*((2*e1/(h(i-1)*h(i)))+(2+x(i)))+1;
    A(i-1,i)=k*(-(2*e1/(h(i)*(h(i-1)+h(i))))+((1+x(i))/(h(i-1)+h(i))));
end
i=2;
   A(i-1,i-1)=k*((2*e1/(h(i-1)*h(i)))+(((1+x(i))+(1+x(i-1)))/(2*h(i-1)))+((2+x(i))+(2+x(i-1)))/2)+(1/2);
    A(i-1,i)=-k*(2*e1/(h(i)*(h(i-1)+h(i))));
 
 i=Nx;   
    A(i-1,i-2)=k*(-(2*e1/(h(i-1)*(h(i-1)+h(i))))-((1+x(i))/(h(i-1)+h(i))));
    A(i-1,i-1)=k*((2*e1/(h(i-1)*h(i)))+(2+x(i)))+1;
     
 for i=3:Nx/2+1
    D(i-1,i-2)=k*(-(2*e2/(h(i-1)*(h(i-1)+h(i))))-(((1+2*x(i))+(1+2*x(i-1)))/(2*h(i-1)))+((2+2*x(i))+(2+2*x(i-1)))/2)+(1/2);
    D(i-1,i-1)=k*((2*e2/(h(i-1)*h(i)))+(((1+2*x(i))+(1+2*x(i-1)))/(2*h(i-1)))+((2+2*x(i))+(2+2*x(i-1)))/2)+(1/2);
    D(i-1,i)=-k*(2*e2/(h(i)*(h(i-1)+h(i))));
end
for i=(Nx/2)+2:Nx-1
    D(i-1,i-2)=k*(-(2*e2/(h(i-1)*(h(i-1)+h(i))))-((1+2*x(i))/(h(i-1)+h(i))));
    D(i-1,i-1)=k*((2*e2/(h(i-1)*h(i)))+(2+2*x(i)))+1;
    D(i-1,i)=k*(-(2*e2/(h(i)*(h(i-1)+h(i))))+((1+2*x(i))/(h(i-1)+h(i))));
end
i=2;
   D(i-1,i-1)=k*((2*e2/(h(i-1)*h(i)))+(((1+2*x(i))+(1+2*x(i-1)))/(2*h(i-1)))+((2+2*x(i))+(2+2*x(i-1)))/2)+(1/2);
    D(i-1,i)=-k*(2*e2/(h(i)*(h(i-1)+h(i))));
 
 i=Nx;   
    D(i-1,i-2)=k*(-(2*e2/(h(i-1)*(h(i-1)+h(i))))-((1+2*x(i))/(h(i-1)+h(i))));
    D(i-1,i-1)=k*((2*e2/(h(i-1)*h(i)))+(2+2*x(i)))+1;
     

  for i=3:Nx/2+1
     B(i-1,i-2)=-k*(x(i)+x(i-1))/4;
     B(i-1,i-1)=-k*(x(i)+x(i-1))/4;
     B(i-1,i)=0;
  end
 for i=(Nx/2)+2:Nx-1
     B(i-1,i-2)=0;
     B(i-1,i-1)=-k*x(i);
     B(i-1,i)=0;
 end
  
 i=2;
     B(i-1,i-1)=-k*(x(i)+x(i-1))/4;
     B(i-1,i)=0;
     
 i=Nx;   
     B(i-1,i-2)=0;
     B(i-1,i-1)=-k*x(i);
    
     M=[A,B;B,D];
  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %for j=2:Nt  
for i=3:(Nx/2)+1
             b1(i-1)=0.5*k*(f1(i-1,j)+f1(i,j))+0.5*(U1(i-1,j-1)+U1(i,j-1));
             b2(i-1)=0.5*k*(f2(i-1,j)+f2(i,j))+0.5*(U2(i-1,j-1)+U2(i,j-1));
end
     for i=(Nx/2)+2:Nx-1
         b1(i-1)=k*f1(i,j)+U1(i,j-1);
         b2(i-1)=k*f2(i,j)+U2(i,j-1);
     end
      b1(1)=0.5*k*(f1(1,j)+f1(2,j))+0.5*(U1(2,j-1));
      b2(1)=0.5*k*(f2(1,j)+f2(2,j))+0.5*(U2(2,j-1));
      b1(Nx-1)=k*f1(Nx+1,j);
      b2(Nx-1)=k*f2(Nx+1,j);
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
     if(l==1)
           F1=U1;
          F2=U2;
    end
       if(l==2)
          G1=U1;
          G2=U2;
       end
       
end
  %%%%%%%%%%%%%%%%%%%%%
%    err1=abs(F1-G1(1:2:Nx+1,1:2:Nt+1));
%     err_sup1=max(err1);
%          max_err1=max(err_sup1);
%       error1(N_counter)=max_err1;
%        %error1(N_counter)=max_err1;
%           
%            %err2=abs(U2-U2_exct);
%         err2=abs(F2-G2(1:2:Nx+1,1:2:Nt+1));
%         err_sup2=max(err2);
%          max_err2=max(err_sup2);
%        error2(N_counter)=max_err2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
 
   erroriginal1=abs(F1-G1(1:2:Nx+1,1:2:Nt+1));
   err11=erroriginal1(1:Nx/2+1);
   err12=erroriginal1(Nx/2+1:Nx+1);
    err_sup11=max(err11);
    err_sup12=max(err12);
         max_err11=max(err_sup11);
         max_err12=max(err_sup12);
      %error1(ep_counter,N_counter)=max_err1;
       error11(N_counter)=max_err11;
       error12(N_counter)=max_err12;  
       
           %err2=abs(U2-U2_exct);
        erroriginal2=abs(F2-G2(1:2:Nx+1,1:2:Nt+1));
        err21=erroriginal2(1:Nx/2+1);
        err22=erroriginal2(Nx/2+1:Nx+1);
       err_sup21=max(err21);
       err_sup22=max(err22);
         max_err21=max(err_sup21);
         max_err22=max(err_sup22);
       %error2(ep_counter,N_counter)=max_err2;
       %error2(N_counter)=max_err2;
         error21(N_counter)=max_err21;
       error22(N_counter)=max_err21;  
end
%end
% %end

% %for j=1:ep_counter
% for i=1:N_counter-1
%          cnv_rt1(i)=log2(error1(i)/error1(i+1));
%          %cnv_rt1(i)=log2(error1(i)/error1(i+1));
% end
% %end
% %for j=1:ep_counter
%  for i=1:N_counter-1
%          cnv_rt2(i)=log2(error2(i)/error2(i+1));
%         %cnv_rt2(i)=log2(error2(i)/error2(i+1));
%  end
%end

% cnv_rt1
% cnv_rt2
% error1
% error2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%for j=1:ep_counter
for i=1:N_counter-1
         %cnv_rt1(j,i)=log2(error1(j,i)/error1(j,i+1));
         cnv_rt11(i)=log2(error11(i)/error11(i+1));
         cnv_rt12(i)=log2(error12(i)/error12(i+1));
end
%end
%for j=1:ep_counter
 for i=1:N_counter-1
         %cnv_rt2(j,i)=log2(error2(j,i)/error2(j,i+1));
         cnv_rt21(i)=log2(error21(i)/error21(i+1));
         cnv_rt22(i)=log2(error22(i)/error22(i+1));
 end
%end
% 
cnv_rt11
cnv_rt12
cnv_rt21
cnv_rt22
error11
error12
error21
error22

% max(cnv_rt1)
% max(cnv_rt2)
% max(error1)
% max(error2)