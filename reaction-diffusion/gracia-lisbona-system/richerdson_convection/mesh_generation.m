function [x,h]=mesh_generation(sig1,sig2,Nx)
for i=1:(Nx/8)
      h(i)=(8*sig1)/Nx;
end

for i=(Nx/8)+1:Nx/4
    h(i)=(8*((sig2)-(sig1)))/Nx;
 end
for i=(Nx/4)+1:3*Nx/4
    h(i)=(2*(1-(2*sig2)))/Nx;
end

for i=(3*Nx/4)+1:(7*Nx/8)
      h(i)=(8*((sig2)-(sig1)))/Nx;
end

for i=(7*Nx/8)+1:Nx
    h(i)=(8*sig1)/Nx;
end
    x(1)=0;
for i=1:Nx

    x(i+1)=x(i)+h(i);
end  
end