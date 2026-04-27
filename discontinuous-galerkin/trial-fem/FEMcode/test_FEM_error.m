function error=test_FEM_error(ne,uh)
% error=test_FEM_error(ne,uh)
% *********************************************
% C.E. Powell 
% -----------------
% Code to compute the energy norm error in the piecewise linear FEM solution
% for a test problem (f=1,p=1,q=0).
%
% Inputs : ne (number of elements)
%          uh (FEM solution)
%
% Outputs: error (energy norm error) 


h=(1/ne); xx=0:h:1; 


% FEM solution at vertices
u1s=uh(1:end-1); 
u2s=uh(2:end);

% quadrature weights and points
weights=h.*[1/6;2/3;1/6];  
x_quad=[xx(1:end-1)',(xx(1:end-1)+h/2)',xx(2:end)'];

Ek2=zeros(ne,1); 


% Quadrature - Simpson's Rule
for i=1:3
    Ek2=Ek2+weights(i).*Ek2_eval(x_quad(:,i),u1s./h,u2s./h);
end

error=sqrt(sum(Ek2));

% % Subroutine to evaluate the integrand at quadrature points
 function Ek2=Ek2_eval(x,u1,u2)
 Ek2=(0.5-x+u1-u2).^2;
  return

