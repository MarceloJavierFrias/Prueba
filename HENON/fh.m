function F=fh(t,x)

	global N p
	%Se pone a cero la funcion con la dimension N adecuada
	
	F = zeros(N,1);

	% Definir los parametros

	ra    = p(1);
	ba    = p(2);
	sigma = p(3);

	%Ecuaciones diferenciales del sistema
	F(1) = 1/(sigma*(x(2)-t));
	F(2) = ((ra-x(3))*t-x(2))/(sigma*(x(2)-t));
	F(3) = (t*x(2)-ba*x(3))/(sigma*(x(2)-t));

end
