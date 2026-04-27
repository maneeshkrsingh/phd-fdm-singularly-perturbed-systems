%function for right hand side........
function yval =sourceff(xval,t,mu)
%yval=(1+exp((xval-1)/epl))*cos(xval)+(1+epl)*(1-exp((xval-1)/epl))*sin(xval);
%yval=1;
%for convection
yval=(-t)*(((-exp(-xval/t))/(t^2*(1-exp(-1/t))))+((-exp(-xval/mu))/(mu^2*(1-exp(-1/mu))))+((pi^2)/2)*sin((pi/2)*xval))-(((exp(-xval/t))/(t*(1-exp(-1/t))))+((exp(-xval/mu))/(mu*(1-exp(-1/mu))))-(pi)*cos((pi/2)*xval))+ 2*(((1-exp(-xval/t))/(1*(1-exp(-1/t))))+((1-exp(-xval/mu))/(1*(1-exp(-1/mu))))-2*sin((pi/2)*xval))-(((1-exp(-xval/mu))/(1*(1-exp(-1/mu))))-xval*exp(xval-1));
%for Reaction
%yval=(-t)*(((1/((t^2)*(1+exp(-1/t))))*(exp(-xval/t)+exp((xval-1)/t)))+((1/((mu^2)*(1+exp(-1/mu))))*(exp(-xval/mu)+exp((xval-1)/mu))))+2*(((1/(1+exp(-1/t)))*(exp(-xval/t)+exp((xval-1)/t)))+((1/(1+exp(-1/mu)))*(exp(-xval/mu)+exp((xval-1)/mu)))+(-2))-(((1/(1+exp(-1/mu)))*(exp(-xval/mu)+exp((xval-1)/mu)))-1);
return
