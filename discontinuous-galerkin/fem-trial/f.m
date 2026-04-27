function f_val = f(e,x)
% y = x;
 

c=1/(1+exp(-1/e));


 f_val=c*(exp(-x/e)+exp(-(1-x)/e))-2;
 
  % y = pi*pi*sin(pi*x);
%y=((pi*pi)*sin(pi*x)+pi*cos(pi*x));
return