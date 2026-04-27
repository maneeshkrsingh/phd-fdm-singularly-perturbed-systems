clc
format short e
m=0.4672;
for  ep_counter=1:2
    e=10^(-ep_counter);
    c1=exp(-1/e);
    c2=1-exp(-1/e);
    
for  N_counter =1:3
        Nx=2^(N_counter+4);
        Nt1=2^(N_counter+4);
        
         sig=min(1/2,m*e*log(Nx1));
 for l=1:2

    %Nx=Nx1*2^(l-1);
    Nt=Nt1*2^(l-1);
            
         
 A=zeros(Nx-1);
b=zeros(Nx-1,1);
D=zeros(Nx-1,1);
U=zeros(Nx+1,Nt+1);
U_etrp=zeros(Nx+1,Nt+1);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Exact solution   
 for j=1:Nt+1
    for i=1:Nx+1
        U_exct(i,j)=exp(-t(j))*(c1+(c2*x(i))-exp((x(i)-1)/e));
    end
 end

%%%%%%%%%%%%construction of finite diffrence scheme
for i=2:Nx-2
        A(i,i-1)=-((2*e*k)/((h(i)+h(i+1))*h(i)))-((k*a(i))/(h(i)));
        A(i,i)=k*((2*e/(h(i)*h(i+1)))+a(i)/h(i))+1;
        A(i,i+1)=-((2*e*k)/((h(i)+h(i+1))*h(i+1))); 
end

i=1;
        A(i,i)=k*((2*e/(h(i)*h(i+1)))+a(i)/h(i))+1;
        A(i,i+1)=-((2*e*k)/((h(i)+h(i+1))*h(i+1)));

i=Nx-1;
        A(i,i-1)=-((2*e*k)/((h(i)+h(i+1))*h(i)))-((k*a(i))/(h(i)));
        A(i,i)=k*((2*e/(h(i)*h(i+1)))+a(i)/h(i))+1;
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
 %%%%%%%%%%%%%%%%%%%%%%%Extrapolation term
U_etrp=2*G(:,1:2:Nt+1)-F;
%%%%%%%%%%%%%%%%%%%%% convergence analysis  
errb=abs(U-U_exct);
      err_supb=max(errb);
        max_errb=max(err_supb);
      errorb(ep_counter,N_counter)=max_errb;
      
    if(N_counter>1)
   cnvb(ep_counter,N_counter-1)=log2(  errorb(ep_counter,N_counter-1)/errorb(ep_counter,N_counter));
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%
   erra=abs(U_etrp-U_exct(:,1:2:Nt+1));
      err_supa=max(erra);
        max_erra=max(err_supa);
      errora(ep_counter,N_counter)=max_erra;
      
    if(N_counter>1)
   cnva(ep_counter,N_counter-1)=log2(  errora(ep_counter,N_counter-1)/errora(ep_counter,N_counter));
   end 
    
end
end
errorb;
errora;

cnvb;

cnva;        
