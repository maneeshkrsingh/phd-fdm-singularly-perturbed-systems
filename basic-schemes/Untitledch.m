clc
clear all
format short e
for  ep_counter=1:3
    e=2^(-2*(ep_counter-1));
 for N_counter =1:2
    Nx=16*2^(N_counter-1);
 sig=0.114
x=zeros(Nx+1,1);
h=zeros(Nx,1);
%%%%%%%%%%%%%%
for i=1:Nx/4
    h(i)=(4*sig)/(Nx);
end
for i=Nx/4+1:3*Nx/4
    h(i)=64*(0.5-2*sig)*((1/Nx^3)+((3*i^2)/Nx^2)+(3/(4*Nx^2))+(3/(16*Nx))-(4.5*i/Nx));
    %h(i)=0.6;
end
for i=3*Nx/4+1:Nx
    h(i)=(4*(sig))/(Nx);
end
% x(1)=0;
%  for i=1:Nx
%     x(i+1)=x(i)+h(i);
%  end
 end
end