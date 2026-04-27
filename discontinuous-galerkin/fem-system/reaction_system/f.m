function y = f(e,x)
   %c=1/(1-exp(-1/e));
   %y = x+c*(1-exp(-(1-x)/e)); %%trial for convection
  %y=c*(1-exp(-(1-x)/e))+(((pi/2)^2)-1)*cos((pi/2)*x)+(pi/2)*sin((pi/2)*x);%%%soucefunction for for convection
 % y=c*(1-exp(-x/e))+c*(1-exp(-(1-x)/e))-exp(x*(x-1))+c*exp(-x/e)+c*exp(-(1-x)/e)+(e^2)*(((2*x-1)^2)+2)*exp(x*(x-1)); %%soucefunction for for reaction
  % y=c/(e)*exp(-x/e)+c/(e)*exp(-(1-x)/e)-(((2*x-1)^2)+2)*exp(x*(x-1));
 
  y=x+e-e;
  %y=e*x/e+x*(x/2+e)-(1/2+e)*c*(exp(-(1-x)/e)-exp(-1/e));
return