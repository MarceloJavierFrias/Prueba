
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Integra las ecuaciones de Lorenz           %%%%
%%%%  Cálcula mapa de Poincaré con el esquema de Henón    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Se borran todas las variables del workspace antes de empezar
clear all
%% Se define dos variables globales: N: orden del sistema
global p  N 

% Orden del sistema
N           = 3;

%% Parametros de la integracion

Rel_Tol     = 1e-10;        % Tolerancia relativa
Abs_Tol     = 1e-12;        % Tolerancia absoluta
Max_Step    = 0.004;         % Paso maximo

%% Todos los parametros se define en la variable option. El Matlab define todos
%% estos parametros por defecto y son todos opcionales para optimizar el
%% calculo.

%options     = odeset('RelTol',Rel_Tol,'AbsTol',Abs_Tol,'Refine',1,'MaxStep ',Max_Step,'Events','events');
options     = odeset('RelTol',Rel_Tol,'AbsTol',Abs_Tol,'Refine',1,'MaxStep',Max_Step);

%% La opcion Events permite controlar la integracion de manera intermedia.
%% Por ejemplo, mi funcion event controla que la variable a_1 no sea mayor
%% que 50, cortando la integracion en caso contrario. La opcion tambien
%% puede ser usada para hacer secciones de Poincare.

%% Directorios de trabajo
%directorio1  = 'D:\Caos_5D_1\Polarizacion izquierda\Integracion ecuaciones';
%directorio2  = 'D:\Caos_5D_1\Polarizacion izquierda\Integracion ecuaciones';
%%  

%% Banderas para seleccionar las opciones del programa: 1-> activa, 0->desactiva
flag_Bifurcacion = 1;      % Calcular un diagrama de bifurcacion (DB)
ver_bifurcacion  = 1;      % Ver el diagrama de bifurcacion del directorio 2

arrancar_archivo = 0;      % Calcular el DB arrancando desde el �ltimo punto del archivo 1
arrancar_manual  = 1;      % Calcular el DB arrancando con cond inicial y parametros dados por el usuario

PI               = 3.14159265359;

%Delta_ra         = 0.1;    % Incremento del parametro nu
T1               = 1000;     % Tiempo de eliminacion de transitorio
T2               = 150;      % Tiempo de integracion
sp               = 0;      % Define el plano de la sección de Poincaré        
  
if arrancar_manual==1
%     disp('Arrancando Manualmente')
%%     Se definen los parametros del sistema

    ra    = 166.1;
    ba     = 8/3;
    sigma  = 10.0; 
            
    p    = [ra  ba  sigma];
       
    X0     = [2.0   -1.0   150.0]; 
    M      = 1;
end
   

% Primer ciclo de integracion

disp(['                                                     '])  
disp(['-----------------------------------------------------']) 
disp(['-----------------------------------------------------']) 
disp(['                                                     ']) 
%disp(['Haciendo Primera integracion del punto ',num2str(j)']) 
disp(['Eliminando Transitorio un tiempo   ', num2str(T1),'  Parametro de control (ra) ',num2str(ra)])
disp(['                                                     ']) 

%  [t,x,TE,YE,IE]  = ode45(@f,[0 T1],X0,options);
[t,x]  = ode45(@f,[0 T1],X0,options);

disp(['Integracion de tiempo  ',num2str(T2),'  Parametro de control (ra) ',num2str(ra)])

%  [t,x,TE,YE,IE]  = ode45(@f,[0 T2],X0,options);
[t,x] = ode45(@f,[0 T2],X0,options);
           
paso = 0

tj = t(:);
xj = x(:,1);   
yj = x(:,2);   
zj = x(:,3);   
                                   
paso = 1
ifi = length(t)
ih = 0;           

for i = 2:ifi
    if xj(i) > sp
       if xj(i-1) < sp
          ih = ih +1;
            Tq = -(xj(i-1) - sp);
            X0 = [10   yj(i-1)   zj(i-1)]; 

            [the,yhe] = ode45(@fh,[0 Tq],X0,options);

            lyhe = length(the);
            y12(ih) = yhe(lyhe,2);
            tc(ih) = tj(i);
            ty(ih,1) = tj(i);
            ty(ih,2) = y12(ih);
       end
    end  
end

   
paso = 2

yp = sort(y12); % Ordena el vector de menor a mayor
tco = sort(tc);

tlen = length(tco)  % Calcula la longitud de tco
tmax = tc(tlen) % Calcula el valor mximo de tc

ipp = length(yp) % Calcula la longitud de yp

ypf = yp(ipp) % Calcula el valor máximo de yp
ypi = yp(1) % Calcula el valor mínimo de yp

mp = zeros(ipp-1,2); % Calcula una matriz de ceros de dimensión ipp-1, 2

for j = 2:ipp
    mp(j,1) = ty(j-1,2);
    mp(j,2) = ty(j,2);
    yx(j) = ty(j-1,2);
    yy(j) = ty(j,2);
end

for jj = 3:ipp
    mp2(jj:1) = y12(jj-2);
    mp2(jj:2) = y12(jj);
    yx2(jj) = y12(jj-2);
    yy2(jj) = y12(jj);
end

for li = 1:800
    ll = 0.1*li;
    rex(li) = li;
    rey(li) = li;
end    
  
save('Poincare_y_r=166,1-Henon.dat', '-ascii', '-double', 'mp');
save('y-vs-t-Poincare_r=166,1-Henon.dat', '-ascii', '-double', 'ty');
filename1 = 'Poincare_y_r=166,1-Henon.pdf';
filename3 = 'y-vs-t-Poincare_r=166,1-Henon.pdf';

figure;
% Establecer propiedades del papel
[wp,hp,xpos,ypos,wbox,hbox,fontsize,fontname,interpr] = papel1x1;

lw = 0.5;    % LineWidth
mk = 1.0;    % MarkerSize
%mk2 = 3.0;    % MarkerSize

xleft1 = ypi; xright1 = ypf;
ydown1 = ypi; yup1 = ypf;

%ttotal = tmax

xleft2 = 0; xright2 = tmax;
ydown2 = ypi; yup2 = ypf;

% Realizar los gráficos
axes('Units','centimeters','Position',[xpos ypos wbox hbox]);







figure(1)
plot(yx,yy,'.r','MarkerSize',lw); hold on
%plot(yx,yy,'-r','LineWidth',lw); hold on
plot(rex,rey,'-k','LineWidth',lw); hold on
 
% Propiedades del gráfico
set(gca,'XLimMode','manual'); set(gca,'YLimMode','manual')
set(gca,'XLim',[xleft1 xright1])
set(gca,'YLim',[ydown1 yup1])

xlabel('$y\left(i\right)$','FontSize',fontsize,'FontName',fontname,'Interpreter',interpr)
ylabel('$y\left(i+1\right)$','FontSize',fontsize,'FontName',fontname,'Interpreter',interpr)
set(gca,'FontSize',fontsize,'FontName',fontname)

%set(gcf, 'renderer', 'painters');

print(gcf, '-dpdf', filename1);

%figure(2)
%plot(yx2,yy2,'.r','MarkerSize',mk2); hold on
%plot(rex,rey,'.k','LineWidth',lw); hold on
 
% Propiedades del gráfico
%set(gca,'XLimMode','manual'); set(gca,'YLimMode','manual')
%set(gca,'XLim',[xleft1 xright1])
%set(gca,'YLim',[ydown1 yup1])

%xlabel('$y\left(i\right)$','FontSize',fontsize,'FontName',fontname,'Interpreter',interpr)
%ylabel('$y\left(i+2\right)$','FontSize',fontsize,'FontName',fontname,'Interpreter',interpr)
%set(gca,'FontSize',fontsize,'FontName',fontname)

%set(gcf, 'renderer', 'painters');

%print(gcf, '-dpdf', filename2);

figure(3)
plot(tc,y12,'.r','LineWidth',mk); hold on
%plot(rex,rey,'.k','LineWidth',lw); hold on
 
% Propiedades del gráfico
set(gca,'XLimMode','manual'); set(gca,'YLimMode','manual')
set(gca,'XLim',[xleft2 xright2])
set(gca,'YLim',[ydown2 yup2])

xlabel('$t$','FontSize',fontsize,'FontName',fontname,'Interpreter',interpr)
ylabel('$y\left(i\right)$','FontSize',fontsize,'FontName',fontname,'Interpreter',interpr)
set(gca,'FontSize',fontsize,'FontName',fontname)

set(gcf, 'renderer', 'painters');

print(gcf, '-dpdf', filename3);


   
  


         
         



