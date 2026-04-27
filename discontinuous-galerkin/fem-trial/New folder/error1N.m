function y = error1N(x1,x2,dU_exct,DU_fem)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
xm = (x1+x2)*0.5;

y=(e^2)*(x2-x1)*((dU_exct(xm)-DU_fem(xm))^2);

end

