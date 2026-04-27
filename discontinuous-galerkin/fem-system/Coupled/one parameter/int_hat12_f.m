function y = int_hat12_f(e,x1,x2)

% This function evaluate \int_{x1}^{x2} f(x) \phi(x) dx,
% where   \phi(x) = (x-x1)/(x2-x1), using the Simpson rule

  xm = (x1+x2)*0.5;
  y = (x2-x1)*(func1(e,x1)*hat2(x1,x1,x2) + 4*func1(e,xm)*hat2(xm,x1,x2)...
               + func1(e,x2)*hat2(x2,x1,x2) )/6;

  return