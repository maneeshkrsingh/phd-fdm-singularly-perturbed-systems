function fval=func_2(x,e1,e2)
c1=1/(1-exp(-1/e1));
c2=1/(1-exp(-1/e2));
fval=-e2*(-(x+2)*exp(x-1)-c2/((e2^2)*exp((x/e2))))-2*(-(x+1)*exp(x-1)-c2/((e2)*exp((x/e2))))-(c1*(1-exp(-x/e1))+c2*(1-exp(-x/e2))-2*sin(pi*x/2))+2*(c2*(1-exp(-x/e2))-x*exp(x-1));