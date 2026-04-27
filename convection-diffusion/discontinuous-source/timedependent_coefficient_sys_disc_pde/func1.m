function f_val=func1(x,t)

if (x<=0.5)
    f_val=-2*(1+x^2)*t^2;
else
    f_val=3*(1+x^2)*t^2;
end
end

