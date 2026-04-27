function yval =sourceff1(xval,t,mu)
%yval=(1+exp((xval-1)/epl))*cos(xval)+(1+epl)*(1-exp((xval-1)/epl))*sin(xval);
%yval=1;
%for convection
yval=(-mu)*(-exp(-xval/mu)/((mu^2)*(1-exp(-1/mu)))-(2*exp(xval-1))-(xval*exp(xval-1)))-2*(exp(-xval/mu)/((mu)*(1-exp(-1/mu)))-(exp(xval-1))-(xval*exp(xval-1)))+4*((1-exp(-xval/mu))/((1-exp(-1/mu)))-(xval*exp(xval-1)))-(((1-exp(-xval/t))/(1*(1-exp(-1/t))))+((1-exp(-xval/mu))/(1*(1-exp(-1/mu))))-2*sin((pi/2)*xval));
%for Reaction
%yval=(-mu)*(((1/((mu^2)*(1+exp(-1/mu))))*(exp(-xval/mu)+exp((xval-1)/mu)))-(((1/(1+exp(-1/t)))*(exp(-xval/t)+exp((xval-1)/t)))+((1/(1+exp(-1/mu)))*(exp(-xval/mu)+exp((xval-1)/mu)))+(-2))+2*(((1/(1+exp(-1/mu)))*(exp(-xval/mu)+exp((xval-1)/mu)))-1));
return