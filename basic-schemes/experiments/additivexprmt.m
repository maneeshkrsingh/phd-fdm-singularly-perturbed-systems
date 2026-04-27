clc
clear
format short e
for ep_counter=1:8
    clear A
    clear B
    clear C
    clear D
    e2=2^(-(ep_counter-1));
    e1=e2*(2^(-2*(ep_counter-1)));
for N_counter =1:4
    Nx1=16*2^(N_counter-1);
    Nt1=(2*4^(N_counter-1));
 Nx1;
 Nt1;
sig2=min(1/4,sqrt(e2)*log(Nx1));
sig1=min(sig2/2,sqrt(e1)*log(Nx1));
% sig2
% sig1
 for l=1:2
     Nx=Nx1*2^(l-1);
    Nt=Nt1*2^(l-1); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U=zeors(2*Nx+2,Nt+1);
U_num=zeros(2*Nx-2,Nt+1);
U1=zeros(Nx+1,Nt+1); %U1(2:Nx,:)=U_num(1:Nx-1,:);
%U_1=zeros(Nx-1,Nt);
U2=zeros(Nx+1,Nt+1); %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:); 
%U_2=zeros(Nx-1,Nt);
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
d=zeros(2*Nx-2,1);
M=zeros(2*Nx-2,2*Nx-2);
x=zeros(Nx+1,1);
h=zeros(Nx,1);
t=zeros(Nt+1,1);
c=zeros(Nx-1,1);
%%%%%%%%%%%%%%%%%%%%%%%
% construction of mesh points along  spatial variable and time variable
 
for i=1:(Nx/8)
      h(i)=(8*sig1)/Nx;
end

for i=(Nx/8)+1:Nx/4
    h(i)=(8*((sig2)-(sig1)))/Nx;
 end
for i=(Nx/4)+1:3*Nx/4
    h(i)=(2*(1-(2*sig2)))/Nx;
end

for i=(3*Nx/4)+1:(7*Nx/8)
      h(i)=(8*((sig2)-(sig1)))/Nx;
end

for i=(7*Nx/8)+1:Nx
    h(i)=(8*sig1)/Nx;
end
h;
 %%%%%%%%%%
for i=1:Nx
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
x;
for j=1:Nt+1
    t(j)=(1*(j-1))/Nt;
end
k=1/Nt;
t;
k;
%%%%%%%%%%%%%%%%%%%
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term
 for j=1:Nt
 for i=1:Nx+1
    U1(i,1)=0;
    U2(i,1)=0;
end

    U1(1,j)=0;
    U1(Nx+1,j)=0;
    
    U2(1,j)=0;
    U2(Nx+1,j)=0;

% for j=1:Nt+1
    for i=1:Nx+1
        f1(i,j)=1;
        f2(i,j)=1;
    end
% end
%%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 % A
 for i=2:Nx-2
     A(i,i-1)=-2*e1*k/((h(i+1)+h(i))*h(i));
     A(i,i)=(2*e1*k/(h(i+1)*h(i)))+k+1;
     A(i,i+1)=(-2*e1*k/((h(i+1)+h(i))*h(i+1)));
 end
 i=1;
      A(i,i)=(2*e1*k/(h(i+1)*h(i)))+k+1;
     A(i,i+1)=(-2*e1*k/((h(i+1)+h(i))*h(i+1)));
     
 i=Nx-1;   
     A(i,i-1)=-2*e1*k/((h(i+1)+h(i))*h(i));
     A(i,i)=(2*e1*k/(h(i+1)*h(i)))+k+1;
%D
 for i=2:Nx-2
     D(i,i-1)=-2*e2*k/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*e2*k/(h(i+1)*h(i)))+k+1;
     D(i,i+1)=(-2*e2*k/((h(i+1)+h(i))*h(i+1)));
 end
 i=1;
      D(i,i)=(2*e2*k/(h(i+1)*h(i)))+k+1;
     D(i,i+1)=(-2*e2*k/((h(i+1)+h(i))*h(i+1)));
     
 i=Nx-1;   
     D(i,i-1)=-2*e2*k/((h(i+1)+h(i))*h(i));
     D(i,i)=(2*e2*k/(h(i+1)*h(i)))+k+1;
     
  for i=1:Nx-1
   c(i,1)=0;
   %c1(i,1)=0;
  end
     B=diag(c);
   C=diag(c);
   M=[A,B;C,D];
   
   
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    for j=2:Nt  
     for i=2:Nx-2
             b1(i)=(k+U1(i,j)+k*U2(i,j));
              b2(i)=(k+U2(i,j)+k*U1(i,j));
              
     end
     
      b1(1)=k;
      b2(1)=k;
      b1(Nx-1)=(k+U1(Nx-1,j)+k*U2(Nx-1,j));
      b2(Nx-1)=(k+U2(Nx-1,j)+k*U1(Nx-1,j));
       
      b=[b1;b2];
         
      d=inv(M)*b;
       %%%%%%%%%%%%%%%%%%%%%%%% numerical solution
%        for i=2:2*Nx
%               U_num(i,j)=d(i-1);    %storing the value of d(each time level) in column of u.
%        end
       U_num(:,j+1)=d;
%        end
   U1(2:Nx,:)=U_num(1:Nx-1,:);
   U2(2:Nx,:)=U_num(Nx:2*Nx-2,:);
    end
    %U1(2:Nx,:)=U_num(1:Nx-1,:);
    %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:);
   U=[U1;U2];
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
 %U_etrp=2*G(:,1:2:Nt+1)-F
  F=[F1,F2];
  G=[G1,G2];
  err1=abs(F1-G1(1:2:Nx+1,1:2:Nt+1));
    err_sup1=max(err1);
         max_err1=max(err_sup1);
       error1(ep_counter,N_counter)=max_err1;
          
          err2=abs(F2-G2(1:2:Nx+1,1:2:Nt+1));
        err_sup2=max(err2);
         max_err2=max(err_sup2);
       error2(ep_counter,N_counter)=max_err2;
          
end
end
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
%     error1
%     error2
%     cnv_rt1
%     cnv_rt2