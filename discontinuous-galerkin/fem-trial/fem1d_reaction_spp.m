
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%         Evaluate the finite element solution                       %
%                                                                    %
%                                                                    %
%             -e^2u''+u=f                                            %
%               u(0)=u(1)=0                                          %
%                                                                    %
%               -e^2u''+2u=f                                         %
%               u(0)=u(1)=0                                          %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all; 
close all; 


for N_counter =1:6
    M=128*2^(N_counter-1);

e=2^(-10);
 c=1/(1+exp(-1/e));
 
 
sig=min(1/4,2.5*e*log(M));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_fem=zeros(M+1,1);
U_exct=zeros(M+1,1);
A =zeros(M-1,M-1); 
F=zeros(M-1,1);  
x=zeros(M+1,1);
h=zeros(M,1);
ERRORON=zeros(M,1);
ERROR1N=zeros(M,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:M/4
    h(i)=4*sig/M;
end



for i=(M/4)+1:3*M/4
    h(i)=2*(1-2*sig)/M;
end

for i=(3*M/4)+1:M
    h(i)=4*sig/M;
end


x(1)=0;
for i=1:M
  x(i+1) = x(i)+h(i);
end

%%%%%%%%%%%%%%%%%%%%mid point %%%%%%%%%%%%%%%
for i=1:M
 xm(i)= (x(i+1)+x(i))/2;
end

% %%%%%%%%%%%%%%%%%%Sparse Matrix %%%%%%%%%
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A(1,1) = 1;
F(1)=0;
A(2,2) = (e^2)/h(1)+2*h(1)/3; 
F(2) = int_hat1_f(e,x(1),x(2));

for i=2:M-2,
  A(i,i) =     A(i,i) + (e^2)/h(i)+2*h(i)/3;
  A(i,i+1)   = A(i,i+1) - (e^2)/h(i)+2*h(i)/6;
  A(i+1,i)   = A(i+1,i) - (e^2)/h(i)+2*h(i)/6;
  A(i+1,i+1) = A(i+1,i+1) + (e^2)/h(i)+2*h(i)/3;
  F(i)       = F(i) + int_hat2_f(e,x(i),x(i+1));
  F(i+1)     = F(i+1) + int_hat1_f(e,x(i),x(i+1));
end

  A(M-1,M-1) = A(M-1,M-1) + (e^2)/h(M-1)+2*h(M-1)/3;
  F(M-1)     = F(M-1) + int_hat2_f(e,x(M-1),x(M));


  U = A\F;
  

%%%%%%%%%%%%% Exact Solution %%%%%%%%%% 
  for i=1:M+1
   %U_exct(i)=x(i)-(exp((x(i)-1)/e)-exp(-(x(i)+1)/e))/(1-exp(-2/e));
   U_exct(i)=c*(exp(-x(i)/e)+exp(-(1-x(i))/e))-1;
  end
  
  %%%%%%%%%%%%% Exact Solution at mid point %%%%%%%%%% 
  
   for i=1:M
   %U_exctm(i)=x(i)-(exp((xm(i)-1)/e)-exp(-(xm(i)+1)/e))/(1-exp(-2/e));
   U_exctm(i)=c*(exp(-xm(i)/e)+exp(-(1-xm(i))/e))-1;
  end
  
  
  %U_fem(2)=U_exct(2);
  
  
 U_fem(1)=0;
 
for i=1:M-1,
  U_fem(i+1) = fem_soln(x,M,U,x(i));  % Compute FEM solution at x2(i)
end

  U_fem(M+1)=0;

  %%%%%%%%%%%%%%5 At mid point %%%%%%%%%%%%%
  
  

  %%%%%%%%%%%%%%%%%%5 derivative of exact solution %%%%%%%%%%%%%%%%%

for i=1:M
   % dU_exct(i)=1-(1/e)*(exp(((x(i))-1)/e)+exp(-((x(i))+1)/e))/(1-exp(-2/e));
   dU_exct(i)=c*(-exp(-x(i)/e)+exp(-(1-x(i))/e))/e;
end
dU_exct(M+1)=0;

%%%%%%%%%%%%%%%%%%%% At mid point %%%%%%%%%%%%%%%%%%%

for i=1:M
    %dU_exctm(i)=1-(1/e)*(exp(((xm(i))-1)/e)+exp(-((xm(i))+1)/e))/(1-exp(-2/e));
    dU_exctm(i)=c*(-exp(-xm(i)/e)+exp(-(1-xm(i))/e))/e;
end
dU_exctm(M+1)=0;

%%%%%%%%%%%%%%%%%%%% difference operator %%%%%%%%%%%%%%%%%
for i=1:M/2
%     DU_fem(i)=2*(((U_fem(i+1)+U_fem(i))/2)-U_fem(i))/h(i);
       % DU_fem(i)=(U_fem(i+1)-U_fem(i))/h(i);
       DU_fem(i)=dU_exctm(i);
end

% DU_fem(1)=dU_exct(1);
for i=M/2:M-1
%     DU_fem(i)=2*(((U_fem(i+1)+U_fem(i))/2)-U_fem(i))/h(i);
       % DU_fem(i)=(U_fem(i+1)-U_fem(i))/h(i);
        DU_fem(i+1)=(U_fem(i+2)-U_fem(i))/(h(i)+h(i+1));
end
DU_fem(M+1)=0;

%%%%%%%%%%%%%%%%%%%%%%% Error Calculation %%%%%%%%%%%%%%%%%%%%
ERRORON(1)=0;
 for i=1:M
 ERRORON(i+1)=ERRORON(i)+(h(i)/6)*(abs((U_exct(i)-U_fem(i)))+4*(abs(U_exctm(i)-(U_fem(i)+U_fem(i+1))/2))+(abs(U_exct(i+1)-U_fem(i+1))));
 end
 
 
 %for i=1:M-1
% ERRORONSUM=sum(ERRORON);
 %end
 
 ERROR1N(1)=0;
 for i=1:M-1
 ERROR1N(i+1)= (ERROR1N(i))+(((e)^2)*(h(i)^2)*((dU_exctm(i)-DU_fem(i))^2));
 end

 
 
%  for i=(M/2)+1:M-2
%  ERROR1N(i+1)= ERROR1N(i)+(e^2)*(h(i))*((dU_exct(i)-DU_fem(i))^2);
%  end
 
 
 ERROR1NFINAL=max(sqrt(ERROR1N));
 
% ERRORFINAL=ERRORONSUM+ERROR1NSUM;
 
ERROR1Nsum(N_counter)=ERROR1NFINAL;
 
 
err1= abs(U_fem-U_exct);

err1_sup=max(err1);
 
   error1(N_counter)=err1_sup;

% 
 error = norm(U_fem-U_exct,inf);

% ERRORON
 
%error
%ERROR1N
%ERROR1Nsum

end

%error
 %ERRORFINAL
 %ERRORONSUM
% ERROR1Nsum


% for i=1:N_counter-1
%          cnv_rt1(i)=log2(error1(i)/error1(i+1));
%         
% end
% 
 for i=1:N_counter-1
          cnv(i)=log2(ERROR1Nsum(i)/ERROR1Nsum(i+1));
         
 end
 ERROR1Nsum
cnv

 %error1
% cnv_rt1

%plot(x,U_fem,'r', x,U_exct,'b'); 
%hold on
plot(x,DU_fem,'r', x,dU_exct,'b'); 
%hold off
%  return