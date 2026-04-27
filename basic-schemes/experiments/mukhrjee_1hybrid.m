clc
clear all
format short e
m=0.42;
for  ep_counter=1:1
    e=10^(-ep_counter);
    c1=exp(-1/e);
    c2=1-exp(-1/e);
    
    for N_counter =1:4
        Nx=2^(N_counter+3);
        Nt=Nx/0.8;
        
        sig=min(1/2,m*e*log(Nx));
        
A=zeros(Nx-1);
b=zeros(Nx-1,1);
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
      h(i)=(2*(1-sig))/(Nx);
end

for i=(Nx/2)+1:Nx
    h(i)=(2*sig)/(Nx);
 end


x(1)=0;
 for i=1:Nx
    x(i+1)=x(i)+h(i);
end

for j=1:Nt+1
    t(j)=(j-1)/Nt;
end
k=1/Nt;

% numerical approximation of initial and boundary value conditions and
% convection term
% nonfomogeneous term
for i=1:Nx+1
    U(i,1)=c1+(c2*x(i))-exp((x(i)-1)/e);
end
for j=1:Nt+1
    U(1,j)=0;
    U(Nx+1,j)=0;
end
for i=1:Nx
    a(i)=1+x(i)-(x(i)^2);
end
for j=1:Nt+1
    for i=1:Nx+1
        f(i,j)=exp(-t(j))*(-c1+(((x(i)^2)/e)-(x(i)/e)+1)*exp((x(i)-1)/e)+c2-((x(i)^2)*c2));
    end
end

%for j=2:Nt+1
    
for i=2:Nx/2
    A(i,i-1)=-k*((2*e/((h(i))*(h(i+1)+h(i))))+((a(i)+a(i-1)))/2*(h(i-1)))+0.5;
   A(i,i)=k*((2*e/(h(i+1)*h(i)))+((a(i)+a(i-1))/2*(h(i))))+0.5;
    A(i,i+1)=-k*(2*e/(h(i+1)*(h(i+1)+h(i))));
end

for i=(Nx/2)+1:Nx-2
     A(i,i-1)=-k*((2*e/(h(i)*(h(i+1)+h(i))))+(a(i-1)/(h(i)+h(i+1))));
    A(i,i)=k*(2*e/(h(i)*h(i+1)))+1;
    A(i,i+1)=-k*((2*e/(h(i)*(h(i)+h(i+1))))-(a(i-1)/(h(i)+h(i+1))));
end
i=1;
   A(i,i)=k*((2*e/(h(i+1)*h(i)))+((a(i)+a(i+1))/2*(h(i))))+0.5;
    A(i,i+1)=-k*(2*e/(h(i+1)*(h(i+1)+h(i))));
 
 i=Nx-1;   
    A(i,i-1)=-k*((2*e/(h(i)*(h(i+1)+h(i))))-(a(i-1)/(h(i)+h(i+1))));
    A(i,i)=k*(2*e/(h(i)*h(i+1)))+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=2:Nt+1
     for i=3:(Nx/2)+1
             b(i-1)=0.5*k*(f(i-1,j)+f(i,j))+(1/2)*(U(i-1,j-1)+U(i,j-1));
     end
     for i=(Nx/2)+2:Nx-1
         b(i-1)=k*(f(i,j))+U(i,j-1);
     end
      b(1)=0.5*k*(f(1,j)+f(2,j))+0.5*(U(1,j-1)+U(2,j-1));
      b(Nx-1)=k*(f(Nx,j)+U(Nx,j-1));
      
       D=A\b;
      %%%%%%%%%%%%%%%%%%%%%%%%% numerical solution
       for i=2:Nx
              U(i,j)=D(i-1);    %storing the value of d(each time level) in column of u.
       end
    end
for j=1:Nt+1
    for i=1:Nx+1
        U_exct(i,j)=exp(-t(j))*(c1+c2*x(i)-exp((x(i)-1)/e));
    end
end
  err=abs(U-U_exct);
      err_sup=max(err);
        max_err=max(err_sup);
      error(ep_counter,N_counter)=max_err;
      
    if(N_counter>1)
   cnv(ep_counter,N_counter-1)=log2(error(ep_counter,N_counter-1)/error(ep_counter,N_counter));
   end 
        
 end
end

error

cnv

 max_error=max(error)
 max_cnv=max(cnv)

