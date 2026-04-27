function fval=func_1(x,e1,e2)
c1=1/(1-exp(-1/e1));
c2=1/(1-exp(-1/e2));
fval=-e1*(0.5*(pi^2)*sin(pi*x/2)-c1*exp(-(x/e1)-2)-c2/((e2^2)*exp((x/e2))))-(-pi*cos(pi*x/2)+c1*exp(-(x/e1)-1)+c2/((e2)*exp((x/e2))))+2*(c1*(1-exp(-x/e1))+c2*(1-exp(-x/e2))-2*sin(pi*x/2))-(c2*(1-exp(-x/e2))-x*exp(x-1));