% Matthew J. Pabin
% Function to plot streamlines, equipotential lines, and pressure contours
% of airflow over a thin airfoil
%
function PlotThinAirfoil(c,alpha,V_inf,p_inf,rho_inf,N)
% Set up grid 
   x_min = -c/2;
   x_max = c*(3/2);
   y_min = -c/2;
   y_max = c/2;
   bound = [x_min x_max y_min y_max];
    
        n_x = 100; % steps in the x direction
        n_y = 100; % steps in the y direction

        [x,y] = meshgrid(linspace(x_min,x_max,n_x),linspace(y_min,y_max,n_y));
% Circulation
dx = c./N;
x_vortex = linspace(dx/2,c-dx,N);
 
g = 2*alpha*V_inf*sqrt((1-(x_vortex/c))./(x_vortex/c)); %strength
G = g.*dx; %circulation

% radius in terms of x so streamlines can be plotted
r =@(x1) sqrt((x-x1).^2 + (y).^2);

% S.L (psi) for uniform flow
psi_uni = V_inf*(y*cos(alpha)-x*sin(alpha));

% V.P (phi) for uniform flo `w
phi_uni = V_inf*(x*cos(alpha)-y*sin(alpha));

% S.L. (psi) for vortex flow
psi_vor = 0;
        for i = 1:N
            psi_vor = psi_vor + (G(i)*log(r(x_vortex(i))))/(2*pi);
        end

% V.P. (phi) for vortex flow
phi_vor = 0;
        for i = 1:N
            theta = atan2(-y,-x+x_vortex(i));
            phi_vor = phi_vor + -(G(i)*theta)/(2*pi);
        end
% Add uniform and vortex flow

 psi = psi_uni + psi_vor;
 phi = phi_uni + phi_vor;

 % Calculate Pressure
 q_inf = (1/2) * rho_inf * V_inf^2; 
Cp = 1-(gradient(phi,(x_max-x_min)/n_x,(y_max-y_min)/n_y)./V_inf).^2;
 p = p_inf + Cp * q_inf;
 
 % Plot streamlines
 figure
contourf(x,y,psi,80) % stream function
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bound)
ylabel('y')
xlabel('x')
title(['Stream Lines for N = ' num2str(N) ' Vorticies']);

% Plot equipotential lines
 figure
contourf(x,y,phi,80) % stream function
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bound)
ylabel('y')
xlabel('x')
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies']);

% Pressure Contours
 figure
contourf(x,y,p,80) % stream function
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bound)
ylabel('y')
xlabel('x')
title(['Pressure Contour Lines for N = ' num2str(N) ' Vorticies']);

end