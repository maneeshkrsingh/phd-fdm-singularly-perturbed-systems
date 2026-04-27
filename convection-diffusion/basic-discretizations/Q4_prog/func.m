function f_val=func(x,e)
% c1=exp(-1/e);
% c2=1-exp(-1/e);
% f_val=exp(-t)*(-c1*(1+exp(1))-c2*(x^2-1+x*exp(1))+(1+(x^2-x)/e+exp(1))*exp((x-1)/e));
if (x<=1/3)
    f_val=-0.5-x;
else
    f_val=9*(x-1)^2;
end

% if (x<=1/2)
%     f_val=-1;
% else
%     f_val=1;
% end