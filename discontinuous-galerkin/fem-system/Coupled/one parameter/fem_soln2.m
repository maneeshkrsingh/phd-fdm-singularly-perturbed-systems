function y = fem_soln2(x,M,U2,xp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
%                  Evaluate the finite element solution              %
%    Input: x, Nodal points                                          %
%           U, The computed coefficients of FEM solution using       %
%              the hat functions                                     %
%           xp,The x coordinate whether the solution will be         %
%              evaluated.                                            %
%    Output: y, The computed FEM solution                            %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M=64;
%M = length(x);
for i=1:M-1,
  if xp >=x(i) && xp <= x(i+1)
     y = hat2(xp,x(i),x(i+1))*U2(i) + hat1(xp,x(i),x(i+1))*U2(i+1);
     return
  end
end