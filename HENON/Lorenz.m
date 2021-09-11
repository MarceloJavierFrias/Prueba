function xp=Lorenz(t,x)
global p
sigma=p(1);
r=p(2);
b=p(3);
x1=x(1);
x2=x(2);
x3=x(3);
xp=[sigma*(x2-x1); (r-x3)*x1-x2; x1*x2-b*x3];
end

