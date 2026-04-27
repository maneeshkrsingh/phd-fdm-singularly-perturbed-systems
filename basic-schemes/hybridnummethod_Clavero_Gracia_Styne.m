clc
clear all
format short e
m=2;
for e_count=1:1
    e=2^(-3*(e_count+1));
 for N_count=1:2
        Nx1=2^(4+N_count);
        Nt1=2^(3+N_count);
        sig=1-(m*e*log(Nx1));
        
        for l=1:2
            Nx=Nx1*2^(l-1);
            Nt=Nt1*2^(l-1);
            
            
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
      h(i)=(2*sig)/(Nx);
end

for i=(Nx/2)+1:Nx
    h(i)=(2*(1-sig))/(Nx);
end 
%   for i=1:Nx
%     H(i)=(h(i)+h(i+1))/2;
%   end
   
  x(1)=0;
 for i=1:Nx
    x(i+1)=x(i)+h(i);
end

for j=1:Nt+1
    t(j)=((j-1))/Nt;
end
k=1/Nt;
t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% coeefficent terms and nonhomo term
for i=1:Nx
    a(i)=(1+((x(i))^2)+(sin(pi*x(i)))/2);
    %a(i-1/2)=((a(i-1)+a(i))/2);
for j=1:Nt
    %b(i,j)=(1+((x(i))^2)+(sin(pi*t(j)))/2);
    f(i,j)=((x(i))^3)*(1-x(i))+(t(j)*(1-t(j))*sin(pi*t(j)));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%numerical approximation of initial and boundary value conditions and
% nonfomogeneous term

for i=1:Nx+1
    U(i,1)=0;
end
for j=1:Nt+1
    U(1,j)=0;
    U(Nx+1,j)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%construction of finite diffrence scheme
for j=2:Nt+1
    
for i=3:Nx/2
    A(i-1,i-2)=-(e/(h(i-1)*(h(i-1)+h(i))/2))-(((a(i-1)+a(i))/2))/h(i-1)+(1/2)*(1+((x(i-1))^2)+(sin(pi*t(j)))/2);
    A(i-1,i-1)=(2*e/(h(i-1)*h(i)))+(((a(i-1)+a(i))/2))/h(i-1)+(1/2)*(1+((x(i))^2)+(sin(pi*t(j)))/2);
    A(i-1,i)=-(e/(h(i)*(h(i-1)+h(i))/2));
end
for i=(Nx/2)+1:Nx-1
    A(i-1,i-2)=-(e/(h(i-1)*(h(i-1)+h(i))/2))-( a(i)/2*(h(i-1)+h(i))/2);
    A(i-1,i-1)=(2*e/(h(i-1)*h(i)))+(1+((x(i))^2)+(sin(pi*t(j)))/2)+1/(k);
    A(i-1,i)=-(e/(h(i)*(h(i-1)+h(i))/2))+( a(i)/2*(h(i-1)+h(i))/2);
end
i=2;
    A(i-1,i-1)=(2*e/(h(i-1)*h(i)))+(((a(i-1)+a(i))/2))/h(i-1)+(1/2)*(1+((x(i))^2)+(sin(pi*t(j))/2))+1/(2*k);
    A(i-1,i)=-(e/(h(i)*(h(i-1)+h(i))/2));
 
 i=Nx;   
    A(i-1,i-2)=-(e/(h(i-1)*(h(i-1)+h(i))/2))-( a(i)/2*(h(i-1)+h(i))/2);
    A(i-1,i-1)=(2*e/(h(i-1)*h(i)))+(1+((x(i))^2)+(sin(pi*t(j))/2))+1/(k);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % for j=2:Nt+1   %check itNt+1
     for i=3:(Nx/2)+1
             d(i-1)=0.5*(f(i-1,j)+f(i,j))+(1/(2*k))*(U(i-1,j-1)+U(i,j-1))+U(i+1,j-1);
     end
     for i=(Nx/2)+2:Nx-1
         d(i-1)=(f(i,j)+((1/k)*U(i,j-1)))+U(i+1,j-1);
     end
      d(1)=0.5*(f(1,j)+f(2,j))+(1/(2*k))*(U(1,j-1)+U(2,j-1))+U(3,j-1);
     d(Nx-1)=(f(Nx+1,j)+((1/(2*k))*U(Nx,j-1)))+U(Nx+1,j-1);
       A ; 
      D=A\d;
      %%%%%%%%%%%%%%%%%%%%%%%%% numerical solution
       for i=2:Nx
              U(i,j)=D(i-1);    %storing the value of d(each time level) in column of u.
       end
 %  end
end
    if(l==1)
          F=U;
      end
       if(l==2)
          G=U;
       end
   end
 err=abs(F-G(1:2:Nx+1,1:2:Nt+1));
   err_sup=max(err);
        max_err=max(err_sup);
      error(e_count,N_count)=max_err
         if(N_count>1)
         cnv_rt(e_count,N_count-1)=log2(error(e_count,N_count-1)/error(e_count,N_count))
         end
 end
 end

    