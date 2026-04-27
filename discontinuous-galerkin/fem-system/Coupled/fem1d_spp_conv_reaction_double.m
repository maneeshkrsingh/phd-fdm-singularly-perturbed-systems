%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%         Evaluate the finite element solution                       %
%                                                                    %
%                                                                    %
%             -e1u1''+u1'+2u1-u2=f1                                  %
%             -e2u1''+2u2'+4u2-u1=f2                                 %
%                                                                    %
%                                                                    %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all; 
close all; 
%format short e
%function U = fem1d(x)    
for N_counter =1:5
    M1=128*2^(N_counter-1);
%M=8192;    
e2=2^(-6);
 e1=2^(-10);
 c1=1/(1-exp(-1/e1));
  c2=1/(1-exp(-1/e2));
sig2=min(1/2,2.5*e2*log(M1));
sig1=min(sig2/2,2.5*(e1)*log(M1));

for l=1:2
     M=M1*2^(l-1);
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U=zeros(2*M,1);
U1_fem=zeros(M+1,1);
U2_fem=zeros(M+1,1);
U1_exct=zeros(M+1,1);
U2_exct=zeros(M+1,1);
U1=zeros(M,1); %U1(2:Nx,1)=U_num(1:Nx-1,1);

U2=zeros(M,1); %U2(2:Nx,1)=U_num(Nx:2*Nx-2,1); 
A = zeros(M-2,M-2);
B = zeros(M-2,M-2);
C = zeros(M-2,M-2);
D = zeros(M-2,M-2);
F1=zeros(M-2,1);  
F2=zeros(M-2,1);
%P=zeros(2*M,1);
%P=zeros(2*M,2*M);
x=zeros(M+1,1);
%x2=zeros(M+1,1);
h=zeros(M,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(M/2)
      h(i)=(2*(1-sig2))/M;
end

for i=(M/2)+1:3*M/4
    h(i)=(4*(sig2-sig1))/M;
 end
for i=(3*M/4)+1:M
    h(i)=(4*sig1)/M;
end


% for i=1:M
%     h(i)=1/M;
% end

 %%%%%%%%%%
for i=1:M
    x(1)=0;
    x(i+1)=x(i)+h(i);
end


%%%%%%%%%%%%%%%%%numerical approximation of  boundary value conditions 
    U1_fem(1)=0;
    U1_fem(M+1)=0;
    
    U2_fem(1)=0;
    U2_fem(M+1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%Sparse Matrix %%%%%%%%%
%%%%%%%%%%%%%%% Initialize
% A(1,1) = 1;
% B(1,1) = 1;
% C(1,1) = 1;
% D(1,1) = 1;
% F1(1)=0;
% F2(1)=0;
% A(M,M) = 1; 
% B(M,M) = 1; 
% C(M,M) = 1; 
% D(M,M) = 1; 
% F1(M)=0;
% F2(M)=0;

%%%%%%%%%%%%%%%%%%%%%%% A and matrix F1 %%%%%%%%%%%%%%%%%%%%%%%%%%%

A(1,1) = e1/h(1)-1/2+2*h(1)/3; 
F1(1) = int_hat11_f(e1,e2,x(1),x(2));

for i=1:M-3,
  A(i,i) = A(i,i) + e1/h(i)-1/2+2*h(i)/3;
  A(i,i+1)   = A(i,i+1) - e1/h(i)+1/2+2*h(i)/6;
  A(i+1,i)   = A(i+1,i) - e1/h(i)-1/2+2*h(i)/6;
  A(i+1,i+1) = A(i+1,i+1) + e1/h(i)+1/2+2*h(i)/3;
  F1(i)       = F1(i) + int_hat12_f(e1,e2,x(i),x(i+1));
  F1(i+1)     = F1(i+1) + int_hat11_f(e1,e2,x(i),x(i+1));
end

  A(M-2,M-2) = A(M-2,M-2) + e1/h(M-2)+1/2+2*h(M-2)/3;
  F1(M-2)     = F1(M-2) + int_hat12_f(e1,e2,x(M-2),x(M-1));


%%%%%%%%%%%%%%%%%%%%%%%%% D and matrix F2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D(1,1) = e2/h(1)-2/2+4*h(1)/3; 
F2(1) = int_hat21_f(e1,e2,x(1),x(2));

for i=1:M-3,
  D(i,i) = D(i,i) + e2/h(i)-2/2+4*h(i)/3;
  D(i,i+1)   = D(i,i+1) - e2/h(i)+2/2+4*h(i)/6;
  D(i+1,i)   = D(i+1,i) - e2/h(i)-2/2+4*h(i)/6;
  D(i+1,i+1) = D(i+1,i+1) + e2/h(i)+2/2+4*h(i)/3;
  F2(i)       = F2(i) + int_hat22_f(e1,e2,x(i),x(i+1));
  F2(i+1)     = F2(i+1) + int_hat21_f(e1,e2,x(i),x(i+1));
end

  D(M-2,M-2) = D(M-2,M-2) + e2/h(M-2)+2/2+4*h(M-2)/3;
  F2(M-2)     = F2(M-2) + int_hat22_f(e1,e2,x(M-2),x(M-1));

  F=[F1;F2];
%%%%%%%%%%%%%%%%%%%%% B and matrix C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B(1,1) = -h(1)/3; 


for i=1:M-3,
  B(i,i) =     B(i,i) -h(i)/3;
  B(i,i+1)   = B(i,i+1) - h(i)/6;
  B(i+1,i)   = B(i+1,i) - h(i)/6;
  B(i+1,i+1) = B(i+1,i+1) -h(i)/3;
 
end
 B(M-2,M-2) = B(M-2,M-2)-h(M-2)/3;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C(1,1) = -h(1)/3; 


for i=1:M-3,
  C(i,i) = C(i,i) -h(i)/3;
  C(i,i+1)   = C(i,i+1) - h(i)/6;
  C(i+1,i)   = C(i+1,i) - h(i)/6;
  C(i+1,i+1) = C(i+1,i+1) -h(i)/3;
 
end
 C(M-2,M-2) = C(M-2,M-2)-h(M-2)/3;
 
%%%%%%%%%%%%  % Full matrix M %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   P=[A,B;C,D];   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



   U =P\F;
U1(1)=0; 
U2(1)=0;

    for i=2:(M-1)
          U1(i)=U(i-1);
    end
   
      for i=M:2*(M)-3
          U2(i+2-M)=U(i-1);
      end 
  
  U1(M)=0;
  U2(M)=0;
%%%%%%%%%%%%% Exact Solution %%%%%%%%%% 
%  for i=1:M
%         
%      U1_exct(i)=c1*(1-exp(-(1-x(i))/e1))+c2*(1-exp(-(1-x(i))/e2))-2*(1-x(i));
%      U2_exct(i)=c2*(1-exp(-(1-x(i))/e2))-(1-x(i))*exp(-x(i));
%   end
% U1_exct(M+1)=0;
% U2_exct(M+1)=0;  

for i=1:M-1,
  U1_fem(i+1) = fem_soln1(x,M,U1,x(i));  % Compute FEM solution at x(i)
end
  

for i=1:M-1,
  U2_fem(i+1) = fem_soln2(x,M,U2,x(i));  % Compute FEM solution at x(i)
end


  if(l==1)
          U11_fem=U1_fem;
          U12_fem=U2_fem;
  end
       if(l==2)
          U21_fem=U1_fem;
          U22_fem=U2_fem;
       end
%        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Error estimation %%%%%%%%%%%%%%%%%%%%%%%%
% err1= abs(U1_fem-U1_exct);
%err2= abs(U2_fem-U2_exct);
err1= abs(U11_fem-U21_fem(1:2:M+1));
 err2= abs(U12_fem-U22_fem(1:2:M+1));
  err1_sup=max(err1);
  err2_sup=max(err2);
   error1(N_counter)=err1_sup;
   error2(N_counter)=err2_sup;
%error = norm(U_fem-U_exct,inf);

end


for i=1:N_counter-1
         cnv_rt1(i)=log2(error1(i)/error1(i+1));
         cnv_rt2(i)=log2(error2(i)/error2(i+1));
end


error1
error2
cnv_rt1
cnv_rt2
% plot(x,U1_fem,'-r', x,U1_exct,'-b') ;
% hold on
% plot(x,U2_fem,'-m', x,U2_exct,'-c') ;
% hold off