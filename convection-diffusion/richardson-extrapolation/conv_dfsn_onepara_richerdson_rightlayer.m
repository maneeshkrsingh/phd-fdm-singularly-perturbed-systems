clc
clear all
format long e

for ep_counter=1:5
    clear A
    clear B
    clear C
    clear D
    e=2^(-2*(ep_counter+2));
   

for N_counter =1:3
    Nx1=16*2^(N_counter-1);
    Nt1=Nx1/2;
 
sig=min(1/2,2*e*log(Nx1));

for l=1:3
     Nx=Nx1*2^(l-1);
    Nt=Nt1*2^(l-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U1=zeros(Nx+1,Nt+1); %U1(2:Nx,:)=U_num(1:Nx-1,:);

U2=zeros(Nx+1,Nt+1); %U2(2:Nx,:)=U_num(Nx:2*Nx-2,:); 

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
      h(i)=2*(1-sig)/Nx;
end


for i=(Nx/2)+1:Nx
    h(i)=(2*(sig))/Nx;
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
    
     U1(i,1)=0;
     U2(i,1)=0;
   
end
 for j=1:Nt
    U1(1,j)=0;
    U1(Nx+1,j)=0;
    
    U2(1,j)=0;
    U2(Nx+1,j)=0;
 end
for j=1:Nt+1
    for i=1:Nx+1
        
     
        f1(i,j)=(t(j)^3)*(x(i)^2)*(1-t(j))*(1-x(i));
        f2(i,j)=(t(j)^3)*(x(i)^2)*(1-x(i))^2;
      
    end
end


%%%%%%%%%%%%%%construction of finite diffrence matrix A B C D M
 % A
 for j=2:Nt+1
 for i=2:Nx-2
     A(i,i-1)=-2*e*k/((h(i+1)+h(i))*h(i))-7*k/h(i);
     A(i,i)=(2*e*k/(h(i+1)*h(i)))+7*k/h(i)+k*(9+x(i))+1;
     A(i,i+1)=(-2*e*k/((h(i+1)+h(i))*h(i+1)));
 end
 i=1;
     A(i,i)=(2*e*k/(h(i+1)*h(i)))+7*k/h(i)+k*(9+x(i))+1;
     A(i,i+1)=(-2*e*k/((h(i+1)+h(i))*h(i+1)));
     
 i=Nx-1;   
     A(i,i-1)=-2*e*k/((h(i+1)+h(i))*h(i))-7*k/h(i);
     A(i,i)=(2*e*k/(h(i+1)*h(i)))+7*k/h(i)+k*(9+x(i))+1;
     
 for i=2:Nx-2
     D(i,i-1)=-2*e*k/((h(i+1)+h(i))*h(i))-7*k/h(i);
     D(i,i)=(2*e*k/(h(i+1)*h(i)))+7*k/h(i)+k*(5+(x(i)))+1;
     D(i,i+1)=(-2*e*k/((h(i+1)+h(i))*h(i+1)));
 end
 i=1;
      D(i,i)=(2*e*k/(h(i+1)*h(i)))+7*k/h(i)+k*(5+(x(i)))+1;
      D(i,i+1)=(-2*e*k/((h(i+1)+h(i))*h(i+1)));
     
 i=Nx-1;   
     D(i,i-1)=-2*e*k/((h(i+1)+h(i))*h(i))-7*k/h(i);
     D(i,i)=(2*e*k/(h(i+1)*h(i)))+7*k/h(i)+k*(5+(x(i)))+1;
     

  for i=2:Nx-2
     B(i,i-1)=0;
     B(i,i)=-k*8;
     B(i,i+1)=0;
 end
 i=1;
     B(i,i)=-k*8;
     B(i,i+1)=0;
     
 i=Nx-1;   
     B(i,i-1)=0;
     B(i,i)=-k*8;
     
   for i=2:Nx-2
     C(i,i-1)=0;
     C(i,i)=-k*4;
     C(i,i+1)=0;
   end
 i=1;
     C(i,i)=-k*4;
     C(i,i+1)=0;
     
 i=Nx-1;   
     C(i,i-1)=0;
     C(i,i)=-k*4; 
     

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
%end
     if(l==1)
           F1=U1;
          F2=U2;
    end
       if(l==2)
          G1=U1;
          G2=U2;
      U1_etrp=2*G1(1:2:Nx+1,1:2:Nt+1)-F1;
      U2_etrp=2*G2(1:2:Nx+1,1:2:Nt+1)-F2;    
          
   %%%%%%%%%%%%%%%%%%%%%%% Error beofre extrapolation term %%%%%%%%%%%%%%%%

  errb1=abs(F1-G1(1:2:Nx+1,1:2:Nt+1));
    err_supb1=max(errb1);
         max_errb1=max(err_supb1);
       errorb1(ep_counter,N_counter)=max_errb1;
       
  
     errb2=abs(F2-G2(1:2:Nx+1,1:2:Nt+1));
    err_supb2=max(errb2);
         max_errb2=max(err_supb2);
       errorb2(ep_counter,N_counter)=max_errb2;    
       end
   
    if(l==3)
        H1=U1;
        H2=U2;
        U11_etrp=2*H1(1:2:Nx+1,1:2:Nt+1)-G1;
        U12_etrp=2*H2(1:2:Nx+1,1:2:Nt+1)-G2;
   
    end
end
   %%%%%%%%%%%%%%%%%%%%%%% Error after extraploation terem %%%%%%%%%%%%%%
      erra1=abs(U1_etrp-U11_etrp(1:2:Nx/2+1,1:2:Nt/2+1));
    err_supa1=max(erra1);
         max_erra1=max(err_supa1);
       errora1(ep_counter,N_counter)=max_erra1;
       
     erra2=abs(U2_etrp-U12_etrp(1:2:Nx/2+1,1:2:Nt/2+1));
    err_supa2=max(erra2);
         max_erra2=max(err_supa2);
       errora2(ep_counter,N_counter)=max_erra2;      
           
           

end
end

for j=1:ep_counter
    for i=1:N_counter-1
         cnv_rtb1(j,i)=log2(errorb1(j,i)/errorb1(j,i+1));
    end
end
for j=1:ep_counter
    for i=1:N_counter-1
         cnv_rtb2(j,i)=log2(errorb2(j,i)/errorb2(j,i+1));
    end
end

for j=1:ep_counter
    for i=1:N_counter-1
         cnv_rta1(j,i)=log2(errora1(j,i)/errora1(j,i+1));
    end
end
for j=1:ep_counter
    for i=1:N_counter-1
         cnv_rta2(j,i)=log2(errora2(j,i)/errora2(j,i+1));
    end
end


%cnv_rtb1
%cnv_rtb2
cnv_rta1
cnv_rta2
%errorb1
%errorb2
errora1
errora2
% max(cnv_rtb1)
% max(cnv_rtb2)
% max(cnv_rta1)
% max(cnv_rta2)
% max(errorb1)
% max(errorb2)
% max(errora1)
% max(errora2)