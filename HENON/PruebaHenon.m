
clear all
global p

%% Parametros de la integracion

Rel_Tol  = 1e-10;        % Tolerancia relativa
Abs_Tol  = 1e-12;        % Tolerancia absoluta
Max_Step = 0.04;        % Paso maximo

options = odeset('RelTol',Rel_Tol,'AbsTol',Abs_Tol,'Refine',1,'MaxStep',Max_Step);

%%     Se definen los parametros del sistema

r = 166.1;
b = 8/3;
sigma = 10.0;

p  = [sigma; r; b];

%% Condiciones inicales

T1 = 1000;
T2 = 550; 
x0 = [2.0; -1.0; 150.0]; 

[t,x]  = ode45(@Lorenz,[0 T1],x0,options);

x0 = [x(length(t),1); x(length(t),2);x(length(t),3)];

[t,x]  = ode45(@Lorenz,[0 T2],x0,options);


%% Se desarma la solución

tj = t(:);
xj = x(:,1);   
yj = x(:,2);   
zj = x(:,3); 


%%
sp = 0;   % Plano de intersección con la hipersuperficie de Poincare
ih = 0; 
for i = 2:length(t)
    if xj(i) > sp     % Busca el cambio de signo
       if xj(i-1) < sp
           ih = ih +1;
           Tq = -(xj(i-1) - sp);
           x0 = [10   yj(i-1)   zj(i-1)]; 
            
           [the,yhe] = ode45(@henon,[0 Tq],x0,options);

           lyhe = length(the);
           y12(ih) = yhe(lyhe,2);
           z12(ih) = yhe(lyhe,3);
           
           tc(ih) = tj(i);
           
           ty(ih,1) = tj(i);
           ty(ih,2) = y12(ih);
           
           tz(ih,1) = tj(i);
           tz(ih,2) = z12(ih);
           
           
           
       end
    end  
end


yp = sort(y12); % Ordena el vector de menor a mayor
tco = sort(tc);

tlen = length(tco);  % Calcula la longitud de tco
tmax = tc(tlen); % Calcula el valor mximo de tc

ipp = length(yp); % Calcula la longitud de yp

ypf = yp(ipp); % Calcula el valor máximo de yp
ypi = yp(1); % Calcula el valor mínimo de yp

mp = zeros(ipp-1,2); % Calcula una matriz de ceros de dimensión ipp-1, 2

for j = 2:ipp
%     mp(j,1) = ty(j-1,2);
%     mp(j,2) = ty(j,2);
    yx(j) = ty(j-1,2);
    yy(j) = ty(j,2);
    
    zx(j) = tz(j-1,2);
    zz(j) = tz(j,2);
    
end

for jj = 3:ipp
    mp2(jj:1) = y12(jj-2);
    mp2(jj:2) = y12(jj);
    yx2(jj) = y12(jj-2);
    yy2(jj) = y12(jj);
end



figure(1)
plot(yx,yy,'.','MarkerSize',7); 
hold on
plot((0:1:round(max(yy))),(0:1:round(max(yy))))
xlabel('$y\left(i\right)$','Interpreter','latex')
ylabel('$y\left(i+1\right)$','Interpreter','latex')



figure(2)
plot(zx,zz,'.','MarkerSize',7); 
hold on
plot((100:1:round(max(zz))),(100:1:round(max(zz))))
xlabel('$z\left(i\right)$','Interpreter','latex')
ylabel('$z\left(i+1\right)$','Interpreter','latex')
axis([100 160 100 160])
