%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%         Evaluate the finite element solution                       %
%                                                                    %
%                                                                    %
%             -eu''+u'+u=c(1-exp(-(1-x)/e))+x                        %
%               u(0)=u(1)=0                                          %
%                                                                    %
%                                                                    %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all; 
close all; 
format short e
%function U = fem1d(x)    
for N_counter =1:2
    M=128*2^(N_counter-1);
%M=8192;    

 e=10^(-5);
 c=1/(1-exp(-1/e));
sig=2.5*(e)*log(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_fem=zeros(M+1,1);
U_exct=zeros(M+1,1);
A = zeros(M,M);
F=zeros(M,1);   
x=zeros(M+1,1);
x2=zeros(M+1,1);
h=zeros(M,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(M/2)
      h(i)=(2*(1-sig))/M;
end

for i=(M/2)+1:M
   h(i)=(2*sig)/M;
end


% for i=1:M
%     h(i)=1/M;
% end

 %%%%%%%%%%
for i=1:M
    x(1)=0;
    x(i+1)=x(i)+h(i);
end
% %%%%%%%%%%%%%%%%%%Sparse Matrix %%%%%%%%%
%%%%%%%%%%%%%%% Initialize
A(1,1) = 1;
F(1)=0;
A(M,M) = 1; 
F(M)=0;


A(2,2) = e/h(1)-1/2; 
F(2) = int_hat1_f(e,x(1),x(2));

for i=2:M-2,
  A(i,i) = A(i,i) + e/h(i)-1/2;
  A(i,i+1)   = A(i,i+1) - e/h(i)+1/2;
  A(i+1,i)   = A(i+1,i) - e/h(i)-1/2;
  A(i+1,i+1) = A(i+1,i+1) + e/h(i)+1/2;
  F(i)       = F(i) + int_hat2_f(e,x(i),x(i+1));
  F(i+1)     = F(i+1) + int_hat1_f(e,x(i),x(i+1));
end

  A(M-1,M-1) = A(M-1,M-1) + e/h(M-1)+1/2;
  F(M-1)     = F(M-1) + int_hat2_f(e,x(M-1),x(M));


% A(2,2) = e/h(1)-1/2+h(1)/3; 
% F(2) = int_hat1_f(e,x(1),x(2));
% 
% for i=2:M-2,
%   A(i,i) = A(i,i) + e/h(i)-1/2+h(i)/3;
%   A(i,i+1)   = A(i,i+1) - e/h(i)+1/2+h(i)/6;
%   A(i+1,i)   = A(i+1,i) - e/h(i)-1/2+h(i)/6;
%   A(i+1,i+1) = A(i+1,i+1) + e/h(i)+1/2+h(i)/3;
%   F(i)       = F(i) + int_hat2_f(e,x(i),x(i+1));
%   F(i+1)     = F(i+1) + int_hat1_f(e,x(i),x(i+1));
% end
% 
%   A(M-1,M-1) = A(M-1,M-1) + e/h(M-1)+1/2+h(M-1)/3;
%   F(M-1)     = F(M-1) + int_hat2_f(e,x(M-1),x(M));


  U = A\F;

  
%%%%%%%%%%%%% Exact Solution %%%%%%%%%% 
  for i=1:M+1
   %U_exct(i)=c*(1-exp(-(1-x(i))/e))+(x(i)-1);
 % U_exct(i)=c*(1-exp(-(1-x(i))/e))-cos((pi/2)*x(i));
 U_exct(i)=x(i)*(x(i)/2+e)-(1/2+e)*c*(exp(-(1-x(i))/e)-exp(-1/e));
  end
  
  

for i=1:M-1,
  U_fem(i+1) = fem_soln(x,M,U,x(i));  % Compute FEM solution at x(i)
end
  
err= abs(U_fem-U_exct);
 err_sup=max(err);
  error(N_counter)=err_sup;
%error = norm(U_fem-U_exct,inf);

end


for i=1:N_counter-1
         cnv_rt(i)=log2(error(i)/error(i+1));
         %cnv_rt1(i)=log2(error1(i)/error1(i+1));
end


error
cnv_rt
plot(x,U_fem,'-r', x,U_exct,'-b') ;
