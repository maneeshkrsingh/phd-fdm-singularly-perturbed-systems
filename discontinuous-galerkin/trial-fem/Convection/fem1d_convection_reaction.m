%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%         Evaluate the finite element solution                       %
%                                                                    %
%                                                                    %
%             -u''+u'+u=(pi^2)*sin(pi*x)+pi*cos(pi*x)+sin(pi*x)      %
%               u(0)=u(1)=0                                          %
%                                                                    %
%                                                                    %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all; 
close all; 
%function U = fem1d(x)        
M=64;                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_fem=zeros(M+1,1);
U_exct=zeros(M+1,1);
A = sparse(M,M);
F=zeros(M,1);   
x=zeros(M+1,1);
x2=zeros(M+1,1);
h=zeros(M,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(1)=0;
for i=1:M
    h(i)=1/M;
end
for i=1:M
  x(i+1) = x(i)+h(i);
end

% %%%%%%%%%%%%%%%%%%Sparse Matrix %%%%%%%%%
 
         % Initialize
A(1,1) = 1; F(1)=0;
A(M,M) = 1; F(M)=0;
A(2,2) = 1/h(1)+h(1)/6-1/2; 
F(2) = int_hat1_f(x(1),x(2));

for i=2:M-2,
  A(i,i) = A(i,i) + 1/h(i)+h(i)/6-1/2;
  A(i,i+1)   = A(i,i+1) - 1/h(i)+h(i)/6+1/2;
  A(i+1,i)   = A(i+1,i) - 1/h(i)+h(i)/6-1/2;
  A(i+1,i+1) = A(i+1,i+1) + 1/h(i)+h(i)/6+1/2;
  F(i)       = F(i) + int_hat2_f(x(i),x(i+1));
  F(i+1)     = F(i+1) + int_hat1_f(x(i),x(i+1));
end

  A(M-1,M-1) = A(M-1,M-1) + 1/h(M-1)+h(M-1)/6+1/2;
  F(M-1)     = F(M-1) + int_hat2_f(x(M-1),x(M));



% for i=2:M-2
%   A(i,i-1)=-1/h(i)+1/2;
%   A(i,i) =  2/h(i)+1/2;
%   A(i,i+1)   = - 1/h(i)+1/2;
%  
% end
% 
%  i=1;
%   A(i,i) =  2/h(i)+1/2;
%   A(i,i+1)   = - 1/h(i)+1/2;
%      
%  i=M-1;   
%   A(i,i-1)=-1/h(i)+1/2;
%   A(i,i) =  2/h(i)+1/2;
% 
%   %%%%%%%%%%%%%%%%%%%%%%% Right hand Side %%%%%%%%%%%%%%
% for i=2:M-3,
%   
%   F(i)       = F(i) + int_hat2_f(x(i),x(i+1));
%   F(i+1)     = F(i+1) + int_hat1_f(x(i),x(i+1));
% end
% F(1) = int_hat1_f(x(1),x(2));
% F(M-1)     = F(M-1) + int_hat2_f(x(M-1),x(M));

% for i=2:M
%   
%   F(i-1)       =  h(i)*((pi*pi)*sin(pi*x(i))+pi*cos(pi*x(i)));
%  % F(i+1)     = F(i+1) + int_hat1_f(x(i),x(i+1));
% end




  U = A\F;

  
%%%%%%%%%%%%% Exact Solution %%%%%%%%%% 
  for i=1:M+1
   U_exct(i)=sin(pi*x(i));
  end
  
  
% U_fem(1)=0;
 %U_fem(M+1)=0;
for i=1:M-1,
  U_fem(i+1) = fem_soln(x,M,U,x(i));  % Compute FEM solution at x2(i)
end
  

error = norm(U_fem-U_exct,inf);
plot(x,U_fem,'-.', x,U_exct) 
%  return