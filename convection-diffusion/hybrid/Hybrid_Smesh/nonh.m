function f_val=nonh(x,e)
c1=exp(-(1-x)/e);
c2=1-exp(-1/e);
f_val=-e*(((-1/e^2)*c1/c2)+((pi^2)/4)*cos(pi*x/2))+(1+x*(1-x))*((-1/e)*c1/c2)+((pi/2)*sin(pi*x/2));
