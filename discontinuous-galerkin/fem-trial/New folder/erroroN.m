function y = erroroN(x1,x2,U_exct,U_fem)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x2=0.1;
xm = (x1+x2)*0.5;
y=((x2-x1)/6)*((U_exct(x1)-U_fem(x1))+4*(U_exct(xm)-U_fem(xm))+(U_exct(x2)-U_fem(x2)));

end

