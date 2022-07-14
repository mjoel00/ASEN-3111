function   PlotThickAirfoil(c,alpha,V_inf,p_inf,rho_inf,N,xx)
% grid 
   x_min = -c/2;
   x_max = c*(3/2);
   y_min = -c/2;
   y_max = c/2;
   bound = [x_min x_max y_min y_max];
    
        n_x = 100;  
        n_y = 100; 
        
        t = xx/100;
        XB_half1 = linspace(c,0,N/2);
        XB_half2 = flip(XB_half1(1:length(XB_half1)-1));
        
        YB_half1 = (t/0.2)*c .* (  0.2969*sqrt(XB_half1/c) - 0.126*(XB_half1/c) - 0.3516*(XB_half1/c).^2 + ...
            0.2843*(XB_half1/c).^3 - 0.1036*(XB_half1/c).^4); 
        YB_half2 = (t/0.2)*c .* (  0.2969*sqrt(XB_half2/c) - 0.126*(XB_half2/c) - 0.3516*(XB_half2/c).^2 + ...
            0.2843*(XB_half2/c).^3 - 0.1036*(XB_half2/c).^4);
       XB = [XB_half1 XB_half2];
       YB = [-YB_half1 YB_half2];
        
      [x,y] = meshgrid(linspace(x_min,x_max,n_x),linspace(y_min,y_max,n_y));
        
[Cl,Cp,G,X,Y] = vortex_panel(XB,YB,V_inf,rad2deg(alpha),0);


% S.L (psi) for uniform flow
psi_uni = V_inf*(y*cos(alpha)-x*sin(alpha));

% V.P (phi) for uniform flow
phi_uni = V_inf*(x*cos(alpha)-y*sin(alpha));

% S.L. (psi) for vortex flow
psi_vor = 0;
        for i = 1:length(X)
            r = sqrt((x-X(i)).^2 + (y-Y(i)).^2);
            psi_vor = psi_vor + ((G(i)*log(r))/(2*pi));
        end

% V.P. (phi) for vortex flow
phi_vor = 0;
        for i = 1:length(X)
            theta = mod(atan2(y-Y(i),x-X(i)),2*pi);
            phi_vor = phi_vor - ((G(i)*theta)/(2*pi));
        end
% Add uniform and vortex flow

 psi = psi_uni + psi_vor;
 phi = phi_uni + phi_vor;

 % Pressure 
  q_inf = (1/2) * rho_inf * V_inf^2; 
Cp = 1-(gradient(phi,(x_max-x_min)/n_x,(y_max-y_min)/n_y)./V_inf).^2;
 p = p_inf + Cp * q_inf;
 
  % Plot streamlines
figure
contourf(x,y,psi,75) % stream function
hold on 
plot(XB,YB,'k','linewidth',2)
hold on
fill(XB,YB,'w')
%hold on
%plot(XB,-YB,'k','linewidth',2)% airfoil
axis(bound)
ylabel('y')
xlabel('x')
title(['Stream Lines for N = ' num2str(N) ' Vorticies']);

% Plot equipotential lines
figure
contourf(x,y,phi,75) % stream function
hold on 
plot(XB,YB,'k','linewidth',2)
hold on
fill(XB,YB,'w')
%hold on
%plot(XB,-YB,'k','linewidth',2)% airfoil
axis(bound)
ylabel('y')
xlabel('x')
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies']);

% Plot Pressure Contours
figure
contourf(x,y,p,80) % stream function
hold on 
plot(XB,YB,'k','linewidth',2) % airfoil
hold on
fill(XB,YB,'w')
axis(bound)
ylabel('y')
xlabel('x')
title(['Pressure Contour Lines for N = ' num2str(N) ' Vorticies']);

