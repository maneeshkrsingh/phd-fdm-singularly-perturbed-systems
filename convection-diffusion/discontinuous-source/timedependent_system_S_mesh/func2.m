function f_val=func2(x,t)

if (x<=0.5)
    f_val=-(1+(x^2)*t^2);
else
    f_val=(1+(x^3)*t);
end
end
