clc
clear
format short e
for ep_counter=1:10
    clear A
    clear B
    clear C
    clear D
    e2=2^(-(ep_counter-1));
    e1=e2*(2^(-(ep_counter-1)));
for N_counter =1:2
    Nx1=16*2^(N_counter-1);
    Nt1=2*2^(N_counter-1);
    
 sig2=min(1/4,sqrt(e2)*log(Nx1));
 sig1=min(sig2/2,sqrt(e1)*log(Nx1));


 for l=1:2
     Nx=Nx1*2^(l-1);
 for m=1:3    
    Nt=Nt1*2^(m-1); 
    
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
% F11=zeros(Nx+1,Nt+1);
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
%error1=zeros(0.5*Nx+1,0.5*Nt+1);
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
 %%%%%%%%%%
for i=1:Nx
    x(1)=0;
    x(i+1)=x(i)+h(i);
end

for j=1:Nt+1
    t(j)=(1*(j-1))/Nt;
end
k=1/Nt;
%%%%%%%%%%%%%%%%%%%
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term
 
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
 for j=2:Nt
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
   c(i,1)=-k;
  end
     B=diag(c);
   C=diag(c);
   M=[A,B;C,D];
   
   
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for j=2:Nt  
      for i=2:Nx-2
             b1(i)=(k+U1(i+1,j-1));
              b2(i)=(k+U2(i+1,j-1));
     end
     
      b1(1)=(k+U1(2,j-1));
      b2(1)=(k+U2(2,j-1));
      b1(Nx-1)=k+U1(Nx,j-1);
      b2(Nx-1)=k+U2(Nx,j-1);
      b=[b1;b2];
      
      d=inv(M)*b;
      
      for i=2:(Nx)
          U1(i,j)=d(i-1);
      end
   
      for i=Nx+1:2*(Nx)-1
          U2(i-Nx+1,j)=d(i-1);
      end
 end
   
    %U1(2:Nx,:)=U_num(1:Nx-1,:);
    %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:);
   U=[U1;U2];
    if(l==1&&m==1)
          F1=U1;
          F2=U2;
      end
       
       if(l==1&&m==2)
           F11=U1;
           F21=U2;
           temp_Nt=Nt;
       end
       
       if(l==2&&m==2)
          G1=U1;
          G2=U2;
       end
       if(l==2&&m==3)
           G11=U1;
           G21=U2;
          % tempr_Nt=Nt;
       end
       
 end
 end
 
 F1_etrp=2*F11(:,1:2:temp_Nt+1)-F1;
  F2_etrp=2*F21(:,1:2:temp_Nt+1)-F2;
 G1_etrp=2*G11(:,1:2:Nt+1)-G1;
 G2_etrp=2*G21(:,1:2:Nt+1)-G2;

   v1=G1_etrp(1:2:Nx+1,:);
 v2=F1_etrp;
  v1_size=size(v1);
 v2_size=size(v2);
 err1=abs(2*v1(:,1:2:v1_size(2))-v2);
  %err1=abs(2*v1(:,1:2:Nt+1)-v2);
    err_sup1=max(err1);
         max_err1=max(err_sup1);
       error1(ep_counter,N_counter)=max_err1;
          if(N_counter>1)
         cnv_rt1(ep_counter,N_counter-1)=log2(error1(ep_counter,N_counter-1)/error1(ep_counter,N_counter));
          end
  w1=G2_etrp(1:2:Nx+1,:);
 w2=F2_etrp;
  w1_size=size(w1);
 w2_size=size(w2);   
          err2=abs(2*w1(:,1:2:w1_size(2))-w2);
          %err2=abs(F2-G2(1:2:Nx+1,1:2:Nt+1));
        err_sup2=max(err2);
         max_err2=max(err_sup2);
       error2(ep_counter,N_counter)=max_err2;
          if(N_counter>1)
         cnv_rt2(ep_counter,N_counter-1)=log2(error2(ep_counter,N_counter-1)/error2(ep_counter,N_counter));
          end
         
end
end
