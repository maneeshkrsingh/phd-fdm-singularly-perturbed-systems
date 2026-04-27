function y = DU_fem(x1,x2,U_fem)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

xm = (x1+x2)*0.5;

y=(U_fem(xm)-U_fem(xm-((x2-x1)/2)))/(x2-x1)/2;

end

