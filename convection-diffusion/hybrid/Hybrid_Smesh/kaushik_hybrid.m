clc
clear all
format short e
m=4.2;
for e_count=1:10
    e=10^(-(e_count));
  for N_count=1:5
        Nx=2^(3+N_count);
       % Nt=2^(2+N_count);
        k=0.8/Nx;
        Nt=1/k;
        sig=m*e*log(Nx);
    m1=exp(-1/e);
    m2=1-exp(-1/e);    
  
 A=zeros(Nx-1);
d=zeros(Nx-1,1);
D=zeros(Nx-1,1);
U=zeros(Nx+1,Nt+1);
U_exct=zeros(Nx+1,Nt+1);
f=zeros(Nx+1,Nt+1);
a=zeros(Nx+1,1);
x=zeros(Nx+1,1);
h=zeros(Nx,1);
t=zeros(Nt+1,1);
% construction of mesh points along time and spatial variable
for i=1:(Nx/2)
      %h(i)=(2*sig)/(Nx);
       h(i)=(2*(1-sig))/(Nx);
end

for i=(Nx/2)+1:Nx
    %h(i)=(2*(1-sig))/(Nx);
    h(i)=(2*sig)/(Nx);
end 
%   for i=1:Nx
%     H(i)=(h(i)+h(i+1))/2;
%   end
   
  x(1)=0;
 for i=1:Nx
    x(i+1)=x(i)+h(i);
 end
x(Nx+1)=1;

 for j=1:Nt+1
    t(j)=(j-1)/Nt;
end
%k=1/Nt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% coeefficent terms and nonhomo term
a(1)=1;
for i=2:Nx
   a(i)=(1+x(i)*(1-x(i)));
 end
  a(Nx+1)=1;

for j=1:Nt+1
    for i=1:Nx+1
     f(i,j)=exp(-t(j))*(-m1+m2*(1-(x(i))^2)+exp(-(1-x(i))/e)*(1-x(i)*(1-x(i))/e));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%numerical approximation of initial and boundary value conditions 
for j=1:Nt+1
    for i=1:Nx+1
     U_exct(i,j)=exp(-t(j))*(m1+(m2*x(i))-exp(-(1-x(i))/e));
    end
end
for i=1:Nx+1
    U(i,1)=(m1+(m2*x(i))-exp(-(1-x(i))/e));
end
for j=1:Nt+1
    U(1,j)=0;
    U(Nx+1,j)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%construction of finite diffrence scheme
for j=2:Nt+1
for i=3:Nx/2+1
    A(i-1,i-2)=k*(-(2*e/(h(i-1)*(h(i-1)+h(i))))-(((1+x(i-1)*(1-x(i-1)))+(1+x(i)*(1-x(i))))/(2*h(i-1))))+(1/2);
    A(i-1,i-1)=k*((2*e/(h(i-1)*h(i)))+(((1+x(i-1)*(1-x(i-1)))+(1+x(i)*(1-x(i))))/(2*h(i-1))))+(1/2);
    A(i-1,i)=-k*(2*e/(h(i)*(h(i-1)+h(i))));
end
for i=(Nx/2)+2:Nx-1
    A(i-1,i-2)=k*(-(2*e/(h(i-1)*(h(i-1)+h(i))))-(a(i)/(h(i-1)+h(i))));
    A(i-1,i-1)=k*((2*e/(h(i-1)*h(i))))+1;
    A(i-1,i)=k*(-(2*e/(h(i)*(h(i-1)+h(i))))+( a(i)/(h(i-1)+h(i))));
end
i=2;
   A(i-1,i-1)=k*((2*e/(h(i-1)*h(i)))+(((1+x(i-1)*(1-x(i-1)))+(1+x(i)*(1-x(i))))/(2*h(i-1))))+(1/2);
    A(i-1,i)=-k*(e/(h(i)*(h(i-1)+h(i))));
 
 i=Nx;   
    A(i-1,i-2)=k*(-(2*e/(h(i-1)*(h(i-1)+h(i))))-(a(i)/(h(i-1)+h(i))));
    A(i-1,i-1)=k*((2*e/(h(i-1)*h(i))))+1;
%%%%%%%%%%%%%%%%%%%% for j=2:Nt+1   %check itNt+1
     for i=3:(Nx/2)+1
             d(i-1)=0.5*k*(f(i-1,j)+f(i,j))+0.5*(U(i-1,j-1)+U(i,j-1));
     end
     for i=(Nx/2)+2:Nx-1
         d(i-1)=k*f(i,j)+U(i,j-1);
     end
      d(1)=0.5*k*(f(1,j)+f(2,j))+0.5*(U(2,j-1));
     d(Nx-1)=k*f(Nx+1,j);
       A ; 
      D=A\d;
      %%%%%%%%%%%%%%%%%%%%%%%%% numerical solution
       for i=2:Nx
              U(i,j)=D(i-1);    %storing the value of d(each time level) in column of u.
       end
 
end

  %%%%%%%%%%%%%%%%%%%%%
 
   err=abs(U-U_exct);
     err_sup=max(err);
         max_err=max(err_sup);
       error(e_count,N_count)=max_err;
  end
end
for j=1:e_count
    for i=1:N_count-1
         cnv_rt(j,i)=log2(error(j,i)/error(j,i+1));
    end
end
cnv_rt
 error