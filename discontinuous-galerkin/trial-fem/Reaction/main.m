
clear all; 
close all;   % Clear every thing so it won't mess up with other
                        % existing variables.

%%%%%% Generate a triangulation

 %x(1)=0; x(2)=0.1; x(3)=0.3; x(4)=0.333; x(5)=0.5; x(6)=0.75;x(7)=1;

M=16;                                    
%M = length(x);
 x=zeros(M+1,1);
x(1)=0;
for i=1:M
    x(i+1)=1/M;
end
%x=0:0.05:1;

%%%%%% Call fem1d function to get the FEM solution at nodal points.

U = fem1d(x);

%%%%%% Compare errors:
 x2=zeros(M+1,1);
x2(1)=0;
for i=1:M
    x2(i+1)=1/M;
end



%x2 = 0:0.01:1; k2 = length(x2);
for i=1:M+1,
  u_exact(i) = soln(x2(i));
  u_fem(i) = fem_soln(x,U,x2(i));  % Compute FEM solution at x2(i)
end

error = norm(u_fem-u_exact,inf)  % Compute the infinity error
   
plot(x2,u_fem,'-.', x2,u_exact)   % Solid: the exact, %dotted: FEM solution
hold; plot(x,U,'o')              % Mark the solution at nodal points.
xlabel('x'); ylabel('u(x) & u_{fem}(x)'); 
title('Solid line: Exact solution, Dotted line: FEM solution')

figure(2); plot(x2,u_fem-u_exact); title('Error plot')   
xlabel('x'); ylabel('u-u_{fem}');  title('Error Plot')