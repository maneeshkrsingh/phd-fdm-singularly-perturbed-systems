function f_val=func(x,t,e)
c1=exp(-1/e);
c2=1-exp(-1/e);
f_val=exp(-t)*(-c1*(1+exp(1))-c2*(x^2-1+x*exp(1))+(1+(x^2-x)/e+exp(1))*exp((x-1)/e));
