function he=henon(t,x)
global p
sigma=p(1);
r=p(2);
b=p(3);
x1=x(1);
x2=x(2);
x3=x(3);
he1 = 1/(sigma*(x(2)-t));
he2 = ((r-x(3))*t-x(2))/(sigma*(x(2)-t));
he3 = (t*x(2)-b*x(3))/(sigma*(x(2)-t));
he = [he1;he2;he3];
end