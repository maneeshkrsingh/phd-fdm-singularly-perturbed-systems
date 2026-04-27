clc
clear all
format short e
m=1;
for  ep_counter=1:1;

    e=2^(-2*(ep_counter-1));
    c1=exp(-1/e);
    c2=1-exp(-1/e);
    
    for N_counter =1:2
    Nx1=16*2^(N_counter-1);
    Nt1=10*2^(N_counter-1);
    
    sig=min(1/2,m*e*log(Nx1));
    
   for l=1:2

    Nx=Nx1*2^(l-1);
%     Nt=Nt1*2^(l-1);
    Nt=Nt1*2^(l-1);
    
A=zeros(Nx-1);
b=zeros(Nx-1,1);
d=zeros(Nx-1,1);
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
    t(j)=(1*(j-1))/Nt;
end
k=1/Nt;
% numerical approximation of initial and boundary value conditions and
% nonfomogeneous term
for i=1:Nx+1
    U(i,1)=c1+(c2*x(i))-exp((x(i)-1)/e);
end
for j=1:Nt+1
    U(1,j)=0;
    U(Nx+1,j)=0;
end

for j=1:Nt+1
    for i=1:Nx+1
        f(i,j)=exp(-t(j))*(-c1+exp((x(i)-1)/e)+c2-((x(i))*c2));
    end
end
%%%%%%%%%%%%construction of finite diffrence scheme
for i=2:Nx-2
        A(i,i-1)=-((2*e*k)/((h(i)+h(i+1))*h(i)))-((k)/(h(i)));
        A(i,i)=1+((2*e*k)/((h(i)+h(i+1))*h(i+1)))+((2*e*k)/((h(i)+h(i+1))*h(i)))+(k/(h(i)));
        A(i,i+1)=-((2*e*k)/((h(i)+h(i+1))*h(i+1))); 
end

i=1;
        A(i,i)=1+((2*e*k)/((h(i)+h(i+1))*h(i+1)))+((2*e*k)/((h(i)+h(i+1))*h(i)))+(k/(h(i)));
        A(i,i+1)=-((2*e*k)/((h(i)+h(i+1))*h(i+1)));   

i=Nx-1;
        A(i,i-1)=-((2*e*k)/((h(i)+h(i+1))*h(i)))-((k)/(h(i)));
        A(i,i)=1+((2*e*k)/((h(i)+h(i+1))*h(i+1)))+((2*e*k)/((h(i)+h(i+1))*h(i)))+(k/(h(i)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=2:Nt+1   %check itNt+1
     for i=2:Nx-2
             b(i)=(k*f(i+1,j))+U(i+1,j-1);
     end
     
      b(1)=(k*f(2,j))+U(2,j-1);
      b(Nx-1)=(k*f(Nx,j))+U(Nx,j-1);
         
      d=A\b;
      %%%%%%%%%%%%%%%%%%%%%%%%% numerical solution
       for i=2:Nx
              U(i,j)=d(i-1);    %storing the value of d(each time level) in column of u.
       end
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
      error(ep_counter,N_counter)=max_err
         if(N_counter>1)
         cnv_rt(ep_counter,N_counter-1)=log2(error(ep_counter,N_counter-1)/error(ep_counter,N_counter))
         end
end
end