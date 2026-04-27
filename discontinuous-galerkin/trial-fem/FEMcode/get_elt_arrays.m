function [Kks,Mks,bks]=get_elt_arrays(h,p,q,f,ne);
% [Kks,Mks,bks]=get_elt_arrays(h,p,q,f,ne);
% *********************************************
% C.E. Powell 
% -----------------
% Subroutine to compute all the element diffusion and reaction 
% FEM matrices and element vectors, for piecewise linear elements
% on uniform meshes with piecewise constant p,q,f.
%
% Inputs : h  (mesh width)
%          p  (diffusion coefficient - a vector of length ne)
%          q  (reaction coefficient - a vector of length ne)
%          f  (rhs term  - a vector of length ne)
%          ne (number of elements)
%
% Outputs: Kks : diffusion element matrices
%          Mks : reaction element matrices
%          bks : element rhs vectors


% Initialise
Kks = zeros(ne,2,2);     
Mks=zeros(ne,2,2); 
bks=zeros(ne,2); 

% diffusion
Kks(:,1,1)=(p./h);       Kks(:,1,2)=-(p./h);
Kks(:,2,1)=-(p./h);      Kks(:,2,2)=(p./h);

% reaction
Mks(:,1,1)=(q.*h./3);    Mks(:,1,2)=(q.*h./6);
Mks(:,2,1)=(q.*h./6);    Mks(:,2,2)=(q.*h./3);

% rhs 
bks(:,1)  = f.*(h./2); 
bks(:,2)  = f.*(h./2);
