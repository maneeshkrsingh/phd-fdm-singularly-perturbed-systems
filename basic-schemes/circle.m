theta=linspace(0,2*pi,100);
r=input('Enter the radius of the circle:')
x=r*cos(theta);
y=r*sin(theta);
plot(x,y)
xlabel('x')
ylabel('y')
title('circle of 2 radius')
print