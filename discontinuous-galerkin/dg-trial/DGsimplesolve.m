clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of the elliptic problem
% -p’’(x) = f(x) in (0,1)
% p(0) = 1
% p(1) = 0
% with primal discontinuous Galerkin methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [Aglobal,rhsglobal,ysol] = DGsimplesolve(nel,ss,penal)
%
% input variables
nel=2;
ss=1;
penal=1;
%sourcef=0;
%locdim=3;
%%%%%sourcefunction
sourcef = @(x) -(2*x-2*(1-2*x)+4*x.*(x-x.^2))*exp(-x.^2);
%%%%%%%%%%%%%%%%%%basisfunction%%%%%%%%%%%%
%h=1/nel;
%phi0=1;
%phi1=(2/h).*(x-(n+1/2).*h);
%phi2=(2/h.^(2)).*(x-(n+1/2).*h).^2;
% output variables
% Aglobal: global stiffness matrix
% rhsglobal: global right-hand side
% ysol: vector of DG unknowns
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%local matrices%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amat = (nel)*[0 0 0;0 4 0;0 0 16/3];
Bmat = (nel)*[penal 1-penal -2+penal; -ss-penal -1+ss+penal 2-ss-penal; 2*ss+penal 1-2*ss-penal -2+2*ss+penal];
Cmat = (nel)*[penal -1+penal -2+penal; ss+penal -1+ss+penal -2+ss+penal; 2*ss+penal -1+2*ss+penal -2+2*ss+penal];
Dmat = (nel)*[-penal -1+penal 2-penal; -ss-penal -1+ss+penal 2-ss-penal; -2*ss-penal -1+2*ss+penal 2-2*ss-penal];
Emat = (nel)*[-penal 1-penal 2-penal; ss+penal -1+ss+penal -2+ss+penal; -2*ss-penal 1-2*ss-penal 2-2*ss-penal];
F0mat =(nel)*[penal 2-penal -4+penal; -2*ss-penal -2+2*ss+penal 4-2*ss-penal; 4*ss+penal 2-4*ss-penal -4+4*ss+penal];
FNmat =(nel)*[penal -2+penal -4+penal; 2*ss+penal -2+2*ss+penal -4+2*ss+penal; 4*ss+penal -2+4*ss+penal -4+4*ss+penal];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
locdim=3;
%%%%%%%%dimension of global matrices%%%%%%
glodim=nel*locdim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aglobal=zeros(glodim, glodim);
rhsglobal = zeros(glodim,1);
%%%%%%%%%%%  gllobal matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%
%first bolck row
for i=1:3
    for j=1:3
        Aglobal(i,j)=Aglobal(i,j)+Amat(i,j)+F0mat(i,j)+Cmat(i,j);
        je=locdim+j;
        Aglobal(i,je)=Aglobal(i,j)+Dmat(i,j);
    end
end
% %intermeadiate block row
% for m=2:(nel-1)
%   for i=1:locdim
%         ie=i+(m-1)*locdim;
%      for j=1:locdim
%             je=j+(m-1)*locdim;
%             Aglobal(ie,je) = Aglobal(ie,je)+Amat(i,j)+Bmat(i,j)+Cmat(i,j);
% je = j+(m-2)*locdim;
% Aglobal(ie,je) = Aglobal(ie,je)+Emat(i,j);
% je = j+(i)*locdim;
% Aglobal(ie,je) = Aglobal(ie,je)+Dmat(i,j);
%       end; 
%    end
% end

        
%%%%%%%%%%%%%%%%%%%%%%%% right hand%%%%%%%%%%%%%%%%%%%        
% Gauss quadrature weights and points %%%%%%%%%%%%%%%%%%%%
wg(1) = 1.0;
wg(2) = 1.0;
sg(1) = -0.577350269189;
sg(2) = 0.577350269189;      
 % compute right-hand side
rhsglobal(1) = nel*penal;
rhsglobal(2) = nel*penal*(-1) - ss*2*nel;
rhsglobal(3) = nel*penal+ss*4*nel;
for i=1:2
rhsglobal(1)= rhsglobal(1)+ wg(i)*sourcef((sg(i)+1)/(2*nel))/(2*nel);
rhsglobal(2)= rhsglobal(2)+ wg(i)*sg(i)*sourcef((sg(i)+1)/(2*nel))/(2*nel);
rhsglobal(3)= rhsglobal(3)+ wg(i)*sg(i)*sg(i)*sourcef((sg(i)+1)/(2*nel))/(2*nel);
end;  

%intermeadiate block row
for m=2:(nel-1)
  for i=1:locdim
        ie=i+(m-1)*locdim;
     for j=1:locdim
            je=j+(m-1)*locdim;
            Aglobal(ie,je) = Aglobal(ie,je)+Amat(i,j)+Bmat(i,j)+Cmat(i,j);
je = j+(m-2)*locdim;
     Aglobal(ie,je) = Aglobal(ie,je)+Emat(i,j);
je = j+(i)*locdim;
    Aglobal(ie,je) = Aglobal(ie,je)+Dmat(i,j);
      end; 
%%%%%%%%%Compute right hand side      
    for ig=1:2
      rhsglobal(ie)= rhsglobal(ie)+wg(ig)*(sg(ig)^(i-1))*sourcef((sg(ig)+2*(m-1)+1.0)/(2*nel))/(2*nel);
    end; 
   end
end
        
% last block row
for i=1:locdim
    ie = i+(nel-1)*locdim;
for j=1:locdim
      je = j+(nel-1)*locdim;
    Aglobal(ie,je) = Aglobal(ie,je)+Amat(i,j)+FNmat(i,j)+Bmat(i,j);
      je = j+(nel-2)*locdim;
    Aglobal(ie,je) = Aglobal(ie,je)+Emat(i,j);
end; 
% compute right-hand side
for ig=1:2
  rhsglobal(ie)= rhsglobal(ie)+wg(ig)*(sg(ig)^(i-1))*sourcef((sg(ig)+2*(nel-1)+1.0)/(2*nel))/(2*nel);
end;
end; 
 %%%%%%%%%%slove linear system%%%%%%%%%%%%%
 alfasolve=Aglobal\rhsglobal;
 
 %%%%%%%%%%%%%%%%end%%%%%%%%%%%%%%
 %%%%%%%%%%%%%display%%%%%%%%%%%%%%%%
 Aglobal
 rhsglobal
 [pol,x]=polation(alfasolve,nel);
 
 sol=@(x) (1-x)*exp(-x^2);
 for i=1:nel+1
     error(i)=sol(x(i))-pol(x(i));
 end
 format short e
 vpa(max(abs(error)),9)
%  plot(x,error)
        
        
        