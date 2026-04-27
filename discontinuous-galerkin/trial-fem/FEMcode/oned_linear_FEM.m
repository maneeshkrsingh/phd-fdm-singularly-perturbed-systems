%function [uh,A,b,K,M]=oned_linear_FEM(ne,p,q,f)
% [uh,A,b,K,M]=oned_linear_FEM(ne,p,q,f)
% *********************************************
% C.E. Powell 
% -----------------
% Code to compute 1d piecewise linear finite element approximation to 
% reaction-diffusion equation on [0,1] with piecewise constant p,q,f, 
% zero  boundary conditions and uniform meshes.
%
% Inputs : ne (no. of elements)
%          p  (diffusion coefficient - a vector of length ne or a single scalar)
%          q  (reaction coefficient - a vector of length ne or a single scalar)
%          f  (rhs term  - a vector of length ne or a single scalar)
%
% Outputs: uh (coefficients of FEM solution)
%          A  (stiffness matrix, A=K+M)
%          b  (rhs of linear system)
%          K  (diffusion matrix)
%          M  (reaction (ie mass) matrix)

ne=16;
p=1;
q=1;
f=1;
% set-up 1d uniform FE mesh
h=(1/ne); 
xx=0:h:1; 
nvtx=length(xx);
J=ne-1; 
elt2vert=[1:J+1;2:(J+2)]';

% initialise global matrices/vectors
K = sparse(nvtx,nvtx); 
M = sparse(nvtx,nvtx); 
b =zeros(nvtx,1);      

% compute element matrices 
[Kks,Mks,bks]=get_elt_arrays(h,p,q,f,ne);

% Assemble element arrays into global arrays
for row_no=1:2
    nrow=elt2vert(:,row_no);
    for col_no=1:2
        ncol=elt2vert(:,col_no);
        K=K+sparse(nrow,ncol,Kks(:,row_no,col_no),nvtx,nvtx);
        M=M+sparse(nrow,ncol,Mks(:,row_no,col_no),nvtx,nvtx);
    end
    b = b+sparse(nrow,1,bks(:,row_no),nvtx,1);
end

% impose homogeneous boundary condition
K([1,end],:)=[]; K(:,[1,end])=[]; 
M([1,end],:)=[]; M(:,[1,end])=[]; 

% Global system 
A=K+M; b(1)=[]; b(end)=[];

% solve linear system for interior degrees of freedom;
u_int=A\b; 


% plot FEM solution
uh=[0;u_int;0]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Error estimate
% 
% % FEM solution at vertices
% u1s=uh(1:end-1); 
% u2s=uh(2:end);
% 
% % quadrature weights and points
% weights=h.*[1/6;2/3;1/6];  
% x_quad=[xx(1:end-1)',(xx(1:end-1)+h/2)',xx(2:end)'];
% 
% Ek2=zeros(ne,1); 
% % Quadrature - Simpson's Rule
% for i=1:3
%     Ek2=Ek2+weights(i).*Ek2_eval(x_quad(:,i),u1s./h,u2s./h);
% end
% 
% error=sqrt(sum(Ek2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(xx,uh,'-o'); axis square; shg;
xlabel('x'); ylabel('u_h(x)');
%title('Piecewise linear FEM approximation')



