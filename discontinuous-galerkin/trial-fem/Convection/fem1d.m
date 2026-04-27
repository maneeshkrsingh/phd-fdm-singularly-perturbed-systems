
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%         Evaluate the finite element solution                       %
%                                                                    %
%                                                                    %
%             -u''+u=(1+(pi^2))sin(pix)                              %
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
A =zeros(M-1,M-1); 
F=zeros(M-1,1);  
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
 
for i=2:M-2
  A(i,i-1)=-1/h(i)+h(i)/6;
  A(i,i) =  2/h(i)+2*h(i)/3;
  A(i,i+1)   = - 1/h(i)+h(i)/6;
 
end

 i=1;
   A(i,i) =  2/h(i)+2*h(i)/3;
  A(i,i+1)   = - 1/h(i)+h(i)/6;
     
 i=M-1;   
   A(i,i-1)=-1/h(i)+h(i)/6;
  A(i,i) =  2/h(i)+2*h(i)/3;

  %%%%%%%%%%%%%%%%%%%%%%% Right hand Side %%%%%%%%%%%%%%
% for i=2:M-3,
%   
%   F(i)       = F(i) + int_hat2_f(x(i),x(i+1));
%   F(i+1)     = F(i+1) + int_hat1_f(x(i),x(i+1));
% end
% F(1) = int_hat1_f(x(1),x(2));
% F(M-1)     = F(M-1) + int_hat2_f(x(M-1),x(M));

for i=2:M
  
  F(i-1)       =  h(i)*(1+pi*pi)*sin(pi*x(i));
 % F(i+1)     = F(i+1) + int_hat1_f(x(i),x(i+1));
end
%F(1) = h*pi*pi*sin(pi*x(2));
%F(M-1)     = h*pi*pi*sin(pi*x(M));



  U = A\F;
  
 
  

  
  x2(1)=0;
for i=1:M
  x2(i+1) = x2(i)+h(i);
end  
  
%%%%%%%%%%%%% Exact Solution %%%%%%%%%% 
  for i=1:M+1
   U_exct(i)=sin(pi*x2(i));
  end
  
  
 U_fem(1)=0;
 U_fem(M+1)=0;
for i=1:M-1,
  U_fem(i+1) = fem_soln(x,U,x2(i));  % Compute FEM solution at x2(i)
end
  

error = norm(U_fem-U_exct,inf)
plot(x2,U_fem,'-.', x2,U_exct) 
%  return