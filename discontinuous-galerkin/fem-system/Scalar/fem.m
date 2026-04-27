clc
clear all
format long

%function C


    
    e=2^(-6);
    % e=1;
   
   
    c=1/(1-exp(-1/e));

  %  Nx=input(' Enter the step-length  ') ;

for N_counter =1:1
    Nx=16*2^(N_counter-1);
   
 

%sig1=min(1/2,2.5*(e)*log(Nx));
sig1=2.5*(e)*log(Nx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U=zeors(2*Nx+2,1);
% U_num=zeros(2*Nx-2,1);


A=zeros(Nx-1);

b1=zeros(Nx-1,1);
f=zeros(Nx+1,1);
g=zeros(Nx+1,1);
U_exct=zeros(Nx+1,1);
C=zeros(Nx-1,1);
F=zeros(Nx-1,1);           
x=zeros(Nx+1,1);
xm=zeros(Nx+1,1);
h=zeros(Nx,1);

%%%%%%%%%%%%%%%%%%%%%%%


% construction of mesh points along  spatial variable and time variable
 
for i=1:(Nx/2)
      h(i)=(2*(1-sig1))/Nx;
end

for i=(Nx/2)+1:Nx
   h(i)=(2*sig1)/Nx;
end


 %%%%%%%%%%
for i=1:Nx
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
xm(i)=0.5*(x(i+1)+x(i));
x;

% numerical approximation of  boundary value conditions 
%     U1(1)=0;
%     U1(Nx+1)=0;
%     
%     U2(1)=0;
%     U2(Nx+1)=0;


%%%%%%%%%%%%%%%%%%%%%%%%%% Exact Solution %%%%%%%%%%%%%%%%%%

  for i=1:Nx+1
        
     U_exct(i)=c*(1-exp(-(1-x(i))/e))+(x(i)-1);
     
  end

 for i=1:Nx+1
        
     f(i)=1+c*(1-exp(-(1-x(i))/e))+(x(i)-1);
     g(i)=1+c*(1-exp(-(1-xm(i))/e))+(xm(i)-1);
  end
%%%%%%%%%%%%%%%%%%%%%% Construction of Block matrices A B C D %%%%%%%%%%%%%%%
%Matrix A

 for i=2:Nx-2
     A(i,i-1)=-e*(1/h(i))+h(i)/6+1/2;
     A(i,i)=e*((1/h(i))+(1/h(i+1)))+(h(i)+h(i+1))/3;
     A(i,i+1)=-e*(1/h(i+1))+h(i+1)/6+1/2;

 end
 i=1;
     A(i,i)=e*((1/h(i))+(1/h(i+1)))+(h(i)+h(i+1))/3;
     A(i,i+1)=-e*(1/h(i+1))+h(i+1)/6+1/2;
     
 i=Nx-1;   
     A(i,i-1)=-e*(1/h(i))+h(i)/6+1/2;
     A(i,i)=e*((1/h(i))+(1/h(i+1)))+(h(i)+h(i+1))/3;
  
%Matrix D  
  
     
%%%%%%%%%%%%%%%%% Right hand Side %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    
     for i=2:Nx-2
             F(i)   =  (x(i+1)-x(i))*(f(i)*hat2(x(i),x(i),x(i+1))+ 4*g(i)*hat2(xm(i),x(i),x(i+1))+ f(i+1)*hat2(x(i+1),x(i),x(i+1)))/3;
             F(i+1)  = (x(i+1)-x(i))*(f(i)*hat1(x(i),x(i),x(i+1))+ 4*g(i)*hat1(xm(i),x(i),x(i+1))+ f(i+1)*hat1(x(i+1),x(i),x(i+1)))/3;
         
             
             
     end
     
           
      F(1)=(x(2)-x(1))*(f(1)*hat1(x(1),x(1),x(2))+ 4*g(1)*hat1(xm(1),x(1),x(2))+ f(2)*hat1(x(2),x(1),x(2)))/3;
      
      F(Nx-1) = (x(Nx)-x(Nx-1))*(f(Nx-1)*hat2(x(Nx-1),x(Nx-1),x(Nx))+ 4*g(Nx-1)*hat2(xm(Nx-1),x(Nx-1),x(Nx))+ f(Nx)*hat2(x(Nx),x(Nx-1),x(Nx)))/3;
     
      
     
   
 %%%%%%%%%%%%%%%%%%%  
 P=A\F;     
      
     for i=2:(Nx)
          C(i)=P(i-1);
      end
   
     
     
    %  C=[C1;C2];
%return
end

% A = sparse(M,M); F=zeros(M,1);           % Initialize
% A(1,1) = 1; F(1)=0;
% A(M,M) = 1; F(M)=0;
% A(2,2) = 1/h(1); F(2) = int_hat1_f(x(1),x(2));
% 
% for i=2:M-2,
%   A(i,i) = A(i,i) + 1/h(i);
%   A(i,i+1)   = A(i,i+1) - 1/h(i);
%   A(i+1,i)   = A(i+1,i) - 1/h(i);
%   A(i+1,i+1) = A(i+1,i+1) + 1/h(i);
%   F(i)       = F(i) + int_hat2_f(x(i),x(i+1));
%   F(i+1)     = F(i+1) + int_hat1_f(x(i),x(i+1));
% end
% 
%   A(M-1,M-1) = A(M-1,M-1) + 1/h(M-1);
%   F(M-1)     = F(M-1) + int_hat2_f(x(M-1),x(M));
% 
%   U = A\F;




%end

%plot(x,C1')
%hold on 
%plot (x,U1_exct')
%hold off