%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%         Evaluate the finite element solution                       %
%                                                                    %
%                                                                    %
%             -e^2u''+u=f                                              %
%               u(0)=u(1)=0                                          %
%                                                                    %
%                                                                    %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all; 
close all; 
%function U = fem1d(x)  
e=10^(-6);

M=128;    

 
 %c=1/(1-exp(-1/e));
  c=1/(1-exp(-1/e));
sig=0.5*(e)*log(M);
%sig=2*(sqrt(e))*log(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_fem=zeros(M+1,1);
U_exct=zeros(M+1,1);
A = zeros(M,M);
F=zeros(M,1);   
x=zeros(M+1,1);
x2=zeros(M+1,1);
h=zeros(M,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:M/4
       h(i)=(4*sig)/M;
   
end

for i=M/4+1:3*M/4
      
    h(i)=(2*(1-2*sig))/M;
end

for i=(3*M/4)+1:M
   h(i)=(4*sig)/M;
end


% for i=1:M
%     h(i)=1/M;
% end

 %%%%%%%%%%
  x(1)=0;
for i=1:M
    x(i+1)=x(i)+h(i);
end
% %%%%%%%%%%%%%%%%%%Sparse Matrix %%%%%%%%%
%%%%%%%%%%%%%%% Initialize
A(1,1) = 1; F(1)=0;
A(M,M) = 1; F(M)=0;
A(2,2) = (e^2)/h(1)+h(1)/3;
F(2) = int_hat1_f(e,x(1),x(2));

for i=2:M-2,
  A(i,i)     = A(i,i) + (e^2)/h(i)+h(i)/3;
  A(i,i+1)   = A(i,i+1) - (e^2)/h(i)+h(i)/6;
  A(i+1,i)   = A(i+1,i) - (e^2)/h(i)+h(i)/6;
  A(i+1,i+1) = A(i+1,i+1) + (e^2)/h(i)+h(i)/3;
  F(i)       = F(i) + int_hat2_f(e,x(i),x(i+1));
  F(i+1)     = F(i+1) + int_hat1_f(e,x(i),x(i+1));
end

  A(M-1,M-1) = A(M-1,M-1) + (e^2)/h(M-1)+h(M-1)/3;
  F(M-1)     = F(M-1) + int_hat2_f(e,x(M-1),x(M));


  U = A\F;

  
%%%%%%%%%%%%% Exact Solution %%%%%%%%%% 
  for i=1:M+1
   U_exct(i)=c*(1-exp(-x(i)/e))+c*(1-exp(-(1-x(i))/e))-exp(x(i)*(x(i)-1));
   % U_exct(i)=x(i)-c*(exp((x(i)-1)/e)+exp(-(1+x(i))/e));
 % U_exct(i)=c*(1-exp(-(1-x(i))/e))-cos((pi/2)*x(i));
  end
  
  

for i=1:M-1,
  U_fem(i+1) = fem_soln(x,M,U,x(i));  % Compute FEM solution at x(i)
end
  
errorpoint= abs(U_fem-U_exct);
error = norm(U_fem-U_exct,inf)
plot(x,U_fem,'-.', x,U_exct) 
