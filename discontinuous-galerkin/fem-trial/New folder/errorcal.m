%%%%%% Error calculation for discrete energy norm%%%%
xmid(i)=(x(i)+x(i+1))/2;
dU_exct=1-(1/e)*(exp((x(i)-1)/e)+exp(-(x(i)+1)/e))/(1-exp(-2/e));
DU_fem=U_fem(xmid(i))-U_fem(xmid(i)-(h(i)/2));
error1N=error1N(i)+(e^2)*(h(i));