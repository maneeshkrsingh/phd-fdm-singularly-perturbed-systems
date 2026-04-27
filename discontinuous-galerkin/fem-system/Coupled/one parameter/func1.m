
%function f_val=func1(e1,e2,x)
function f_val=func1(e,x)
c=1/(1-exp(-1/e));
%c2=1/(1-exp(-1/e2));

%f_val=exp(x)+e1/e1-e2/e2; %%%%%%%%%%%%double mesh priniciple 
f_val=2+3*c*(1-exp(-(1-x)/e))-8*(1-x)+(1-x)*exp(x); %%%%%%%%%%%%one parameteer
%%%%%%%%%%%%f1=-e1u''_(1)+u'_(1)+2u_(1)-u_(2)%%%%%%%%%%%%%%
 %f_val=(e1/(e2)^2)*c2*(exp(-(1-x)/e2))-(1/e2)*c2*(exp(-(1-x)/e2))+2+2*c1*(1-exp(-(1-x)/e1))+c2*(1-exp(-(1-x)/e2))-4*(1-x)+(1-x)*exp(x);
  
end

