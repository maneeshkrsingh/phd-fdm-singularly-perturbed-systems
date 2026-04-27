function f_val=func22(x,y,t,e)
m1=-exp(-1/e);
m2=-1-m1;
%T1=(1+t-exp(-t));
T2=(1+t^(2)-exp(-t));
Ax=(m1+m2*x+exp(-(1-x)/e));
%Ay=(m1+m2*y+exp(-(1-y)/e));

f_val=-e*T2*Ax*((1/e^2)*exp(-1/e))+T2*Ax*(m2+(1/e)*exp(-1/e))+y*(e*T2*Ax*(1/e^2)+T2*Ax*m2+e*T2*Ax*((1/e^2)*exp(-1/e))-T2*Ax*(m2+(1/e)*exp(-1/e))); 
end