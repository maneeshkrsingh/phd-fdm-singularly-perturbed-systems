function y = dU_exct(x1,x2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
xm = (x1+x2)*0.5;
y=1-(1/e)*(exp(((xm)-1)/e)+exp(-((xm)+1)/e))/(1-exp(-2/e));

end

