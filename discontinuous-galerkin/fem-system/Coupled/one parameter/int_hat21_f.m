function y = int_hat21_f(e,x1,x2)

% This function evaluate \int_{x1}^{x2} f(x) \phi(x) dx,
% where   \phi(x) = (x-x1)/(x2-x1), using the Simpson rule
 
  xm = (x1+x2)*0.5;
  y = (x2-x1)*(func2(e,x1)*hat1(x1,x1,x2) + 4*func2(e,xm)*hat1(xm,x1,x2)...
               + func2(e,x2)*hat1(x2,x1,x2) )/6;

  return
