clc
clear all
format long

%function C


    e2=2^(-2);
    e1=2^(-4);
   
    c1=1/(1-exp(-1/e1));
    c2=1/(1-exp(-1/e2));

  %  Nx=input(' Enter the step-length  ') ;

for N_counter =1:1
    Nx=16*2^(N_counter-1);
   
 
sig2=min(1/2,2.5*e2*log(Nx));
sig1=min(sig2/2,2.5*(e1)*log(Nx));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U=zeors(2*Nx+2,1);
% U_num=zeros(2*Nx-2,1);
U1=zeros(Nx+1,1); %U1(2:Nx,1)=U_num(1:Nx-1,1);

U2=zeros(Nx+1,1); %U2(2:Nx,1)=U_num(Nx:2*Nx-2,1); 

U1_exct=zeros(Nx+1,1);
U2_exct=zeros(Nx+1,1);

A=zeros(Nx-1);
B=zeros(Nx-1);
C=zeros(Nx-1);
D=zeros(Nx-1);
b1=zeros(Nx-1,1);
b2=zeros(Nx-1,1);
P=zeros(2*Nx-2,1);
M=zeros(2*Nx-2,2*Nx-2);
C1=zeros(Nx-1,1);
C2=zeros(Nx-1,1);
x=zeros(Nx+1,1);
h=zeros(Nx,1);

%%%%%%%%%%%%%%%%%%%%%%%


% construction of mesh points along  spatial variable and time variable
 
for i=1:(Nx/2)
      h(i)=(2*(1-sig2))/Nx;
end

for i=(Nx/2)+1:3*Nx/4
    h(i)=(4*((sig2)-(sig1)))/Nx;
 end
for i=(3*Nx/4)+1:Nx
    h(i)=(4*sig1)/Nx;
end


 %%%%%%%%%%
for i=1:Nx
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
x;

% numerical approximation of  boundary value conditions 
%     U1(1)=0;
%     U1(Nx+1)=0;
%     
%     U2(1)=0;
%     U2(Nx+1)=0;


%%%%%%%%%%%%%%%%%%%%%%%%%% Exact Solution %%%%%%%%%%%%%%%%%%

  for i=1:Nx+1
        
     U1_exct(i)=c1*(1-exp(-(1-x(i))/e1))+c2*(1-exp(-(1-x(i))/e2))-2*(1-x(i));
     U2_exct(i)=c2*(1-exp(-(1-x(i))/e2))-(1-x(i))*exp(-x(i));
  end


%%%%%%%%%%%%%%%%%%%%%% Construction of Block matrices A B C D %%%%%%%%%%%%%%%
%Matrix A

 for i=2:Nx-2
     A(i,i-1)=-e1*(1/h(i))+h(i)/3+1/2;
     A(i,i)=e1*((1/h(i))+(1/h(i+1)))+2*(h(i)+h(i+1))/3;
     A(i,i+1)=-e1*(1/h(i+1))+h(i+1)/3+1/2;

 end
 i=1;
     A(i,i)=e1*((1/h(i))+(1/h(i+1)))+2*(h(i)+h(i+1))/3;
     A(i,i+1)=-e1*(1/h(i+1))+h(i+1)/3+1/2;
     
 i=Nx-1;   
     A(i,i-1)=-e1*(1/h(i))+h(i)/3+1/2;
     A(i,i)=e1*((1/h(i))+(1/h(i+1)))+2*(h(i)+h(i+1))/3;
  
%Matrix D  
  
  for i=2:Nx-2
     D(i,i-1)=-e2*(1/h(i))+2*h(i)/3+1/2;
     D(i,i)=e2*((1/h(i))+(1/h(i+1)))+4*(h(i)+h(i+1))/3;
     D(i,i+1)=-e2*(1/h(i+1))+2*h(i+1)/3+1/2;

 end
 i=1;
     D(i,i)=e2*((1/h(i))+(1/h(i+1)))+4*(h(i)+h(i+1))/3;
     D(i,i+1)=-e2*(1/h(i+1))+2*h(i+1)/3+1/2;
     
 i=Nx-1;   
     D(i,i-1)=-e2*(1/h(i))+2*h(i)/3+1/2;
     D(i,i)=e2*((1/h(i))+(1/h(i+1)))+4*(h(i)+h(i+1))/3;
     
    
 % Matrix B
  
 for i=2:Nx-2
     B(i,i-1)=-h(i)/6;
     B(i,i)=-(h(i)+h(i+1))/3;
     B(i,i+1)=-h(i+1)/6;
 end
 i=1;
     B(i,i)=-(h(i)+h(i+1))/3;
     B(i,i+1)=-h(i+1)/6;
     
 i=Nx-1;   
     B(i,i-1)=-h(i)/6;
     B(i,i)=-(h(i)+h(i+1))/3;
     
     
 % Matrix C
  
 for i=2:Nx-2
     C(i,i-1)=-h(i)/6;
     C(i,i)=-(h(i)+h(i+1))/3;
     C(i,i+1)=-h(i+1)/6;
 end
 i=1;
     C(i,i)=-(h(i)+h(i+1))/3;
     C(i,i+1)=-h(i+1)/6;
     
 i=Nx-1;   
     C(i,i-1)=-h(i)/6;
     C(i,i)=-(h(i)+h(i+1))/3;
     
 % Full matrix M
   M=[A,B;C,D];   
     
%%%%%%%%%%%%%%%%% Right hand Side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    
     for i=2:Nx-2
             b1(i)=(1/2)*(h(i)+h(i+1))*func1(e1,e2,x(i));
             b2(i)=(1/2)*(h(i)+h(i+1))*func2(e1,e2,x(i));
     end
     
      b1(1)=(1/2)*(h(1)+h(2))*func1(e1,e2,x(1));
      b2(1)=(1/2)*(h(1)+h(2))*func2(e1,e2,x(1));
      b1(Nx-1)=(1/2)*(h(Nx-1)+h(Nx))*func1(e1,e2,x(Nx-1));
      b2(Nx-1)=(1/2)*(h(Nx-1)+h(Nx))*func1(e1,e2,x(Nx-1));
      
      b=[b1;b2];
   
 %%%%%%%%%%%%%%%%%%%  
 P=M\b;     
      
     for i=2:(Nx)
          C1(i)=P(i-1);
      end
   
      for i=Nx+1:2*(Nx)-1
          C2(i+1-Nx)=P(i-1);
      end 
     
    %  C=[C1;C2];
%return
end

%end

%plot(x,C1')
%hold on 
%plot (x,U1_exct')
%hold off

