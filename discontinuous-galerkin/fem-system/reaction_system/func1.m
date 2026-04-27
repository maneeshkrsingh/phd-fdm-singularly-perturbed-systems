
function f_val=func1(e1,e2,x)

c1=1/(1+exp(-1/e1));
c2=1/(1+exp(-1/e2));


%f_val=-c2*(exp(-x/e2)+exp(-(1-x)/e2))+c1*(exp(-x/e1)+exp(-(1-x)/e1))-2;

f_val=-c1*(exp(-x/e1)+exp(-(1-x)/e1))-((e1)^2/(e2)^2)*c2*(exp(-x/e2)+exp(-(1-x)/e2))+2*c1*(exp(-x/e1)+exp(-(1-x)/e1))+2*c2*(exp(-x/e2)+exp(-(1-x)/e2))-4-c2*(exp(-x/e2)+exp(-(1-x)/e2))+1;
  
end

