function f_val=func2(e1,e2,x)
%function f_val=func2(e,x)
%c=1/(1-exp(-1/e));
c1=1/(1-exp(-1/e1));
c2=1/(1-exp(-1/e2));

% fe1(x)=exp(-(1-x)/e1);
% fe2(x)=exp(-(1-x)/e2);

%f_val=cos(x)+e1/e1-e2/e2; %%%%%%%%%%%%%%% double mesh principle 
% f_val=-e2*(x+1)*exp(x)-(1/e)*c*(exp(-(1-x)/e))+2*x*exp(x)+2*c*(1-exp(-(1-x)/e))-4*(1-x)*exp(x)+2*(1-x); %%%%%%%%%% one parameter
%%%%%%%%%%%%f2=-e2u''_(2)+2u'_(2)+4u_(2)-u_(1)%%%%%%%%%%%%%%
    %f_val=(-e2*(x+1)*exp(x)-(1/e2)*c2*(exp(-(1-x)/e2))+2*x*exp(x)+3*c2*(1-exp(-(1-x)/e2))-4*(1-x)*exp(x)-c1*(1-exp(-(1-x)/e1))+2*(1-x));
    
    f_val=-(1/e2)*c2*(exp(-(1-x)/e2))+2+3*c2*(1-exp(-(1-x)/e2))-4*(1-x)-c1*(1-exp(-(1-x)/e1))+2*(1-x);
    
    
  %  f_val=(-e2*(pi^2)*sin(pi*(1-x)/2)+(e2/e1^2)*c1*(exp(-(1-x)/e2))-(2/e1)*c1*(exp(-(1-x)/e1))-(1/e2)*c2*(exp(-(1-x)/e2))+4*pi*cos(pi*(1-x)/2)+3*c1*(1-exp(-(1-x)/e1))+3*c2*(1-exp(-(1-x)/e2))-8*sin(pi*(1-x)/2)+2*(1-x));

end
