clc
clear all
format short e

    e2=2^(-5);
    e1=2^(-7);
  
    c1=1/(1-exp(-1/e1));
    c2=1/(1-exp(-1/e2));

%for N_counter =1:2
    Nx=64;
  
    
 
sig2=min(1/2,2*e2*log(Nx));
sig1=min(sig2/2,2*(e1)*log(Nx));



% construction of mesh points 
 
for i=1:(Nx/2)
      h(i)=2*(1-sig2)/Nx;
end

for i=(Nx/2)+1:3*Nx/4
    h(i)=(4*((sig2)-(sig1)))/Nx;
 end
for i=(3*Nx/4)+1:Nx
    h(i)=(4*(sig1))/Nx;
end


 %%%%%%%%%%
for i=1:Nx-1
    x(1)=0;
    x(i+1)=x(i)+h(i);
    x(Nx+1)=1;
end
x;


 for i=1:Nx+1

     U1(i)=(c1*(1-exp(-(1-x(i))/e1))-(1-x(i)));
      U2(i)=(c2*(1-exp(-(1-x(i))/e2))-(1-x(i)));
  
end

figure
subplot(2,2,1)
 plot(x,U1)
 hold on
 plot(x,U2,'r')
 hold off
%plot(x,y1)
title('')

subplot(2,2,2)
 plot(x,U1)
 hold on
 plot(x,U2,'r')
 hold off
title('')

subplot(2,2,3)
 plot(x,U1)
 hold on
 plot(x,U2,'r')
 hold off
title('')


