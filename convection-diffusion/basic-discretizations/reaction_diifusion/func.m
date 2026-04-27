function f_val=func(x,e)
% c1=exp(-1/e);
% c2=1-exp(-1/e);
% f_val=exp(-t)*(-c1*(1+exp(1))-c2*(x^2-1+x*exp(1))+(1+(x^2-x)/e+exp(1))*exp((x-1)/e));
% % if (x<=0.25)
% %     f_val=-4*x;
% % end
% % if (0.25<x<=0.4)
% %     f_val=-1;
% % end
% % if (0.4<x<=0.5)
% %     f_val=1;
% % end
% % if (0.5<x<=1)
% %     f_val=2-2*x;
% % end
if (x<=0.5)
    f_val=2*(1+x^2);
else
    f_val=3*(1+x^2);
end
end

% if (x<=1/2)
%     f_val=-1;
% else
%     f_val=1;
% end