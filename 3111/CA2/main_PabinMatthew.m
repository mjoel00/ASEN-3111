% Matthew J. Pabin
% ASEN 3111: Aerodynamics
% Computational Assignment #2
% 2/7/2022

clear
clc
close all

%% Known Values
c = 2;   %m
alpha = deg2rad(9);  %rad
V_inf  = 60;    %m/s
rho_inf = 1;   %kg/m^3
p_inf = 85.5 * 10^3;  %Pa
N = 500; % # of vorticies
xx = 12;

%% Call plotthin & conduct study
fprintf('------ Problem 1 -----')
% Call function to plot thin airfoil
PlotThinAirfoil(c,alpha,V_inf,p_inf,rho_inf,N);

% Conduct a study
N10 = 10;
PlotThinAirfoil(c,alpha,V_inf,p_inf,rho_inf,N10);
N100 = 100;
PlotThinAirfoil(c,alpha,V_inf,p_inf,rho_inf,N100);
N1000 = 1000;
PlotThinAirfoil(c,alpha,V_inf,p_inf,rho_inf,N1000);
fprintf('\n')
fprintf('\n')
fprintf('Figures 1-3 show the flow and pressure contour for my base value N=500 vorticies')
fprintf('\n')
fprintf('\n')
fprintf('Conduct a study: Figures 4-12 show the flow and pressure contour around a THIN airfoil when N is varied at 10,100, and 1000')
fprintf('\n')
fprintf('\n')
fprintf('\n')
%% Call plotthick & conduct study
fprintf('------ Problem 2 -----')
% Call function to plot thick airfoil
PlotThickAirfoil(c,alpha,V_inf,p_inf,rho_inf,N,xx)

% Conduct a study
PlotThickAirfoil(c,alpha,V_inf,p_inf,rho_inf,N10,xx)
PlotThickAirfoil(c,alpha,V_inf,p_inf,rho_inf,N100,xx)
PlotThickAirfoil(c,alpha,V_inf,p_inf,rho_inf,N1000,xx)
fprintf('\n')
fprintf('\n')
fprintf('Figures 13-15 show the flow and pressure contour for my base value N=500 vorticies')
fprintf('\n')
fprintf('\n')
fprintf('Conduct a Study: Figures 16-24 show the flow and pressure contour around a THICK airfoil when N is varied at 10,100, and 1000')
fprintf('\n')
fprintf('\n')
fprintf('\n')

%% Convergence of Cl vs Alpha
Cl = zeros(1,1000);
for  N = 10:500
        t = xx/100;
    
        XB_half1 = linspace(c,0,N/2);
        XB_half2 = flip(XB_half1(1:length(XB_half1)-1));
        
        YB_half1 = (t/0.2)*c .* (  0.2969*sqrt(XB_half1/c) - 0.126*(XB_half1/c) - 0.3516*(XB_half1/c).^2 + ...
            0.2843*(XB_half1/c).^3 - 0.1036*(XB_half1/c).^4); 
        YB_half2 = (t/0.2)*c .* (  0.2969*sqrt(XB_half2/c) - 0.126*(XB_half2/c) - 0.3516*(XB_half2/c).^2 + ...
            0.2843*(XB_half2/c).^3 - 0.1036*(XB_half2/c).^4);
       XB = [XB_half1 XB_half2];
       YB = [-YB_half1 YB_half2];
       
      [CL(N),~,~,~,~] = vortex_panel(XB,YB,V_inf,alpha,0); 
      
end
 CL = CL(10:end);
 CL_a = CL(end);
 N = 10:500;
 N_max = 500;
 plot(N,CL,'r','LineWidth',2)
 xlabel('# of Discrete Vorticies')
 ylabel('Sectional Lift Coefficient')
 title('Convergence of Cl')
 fprintf('------ Problem 3 ------')
 fprintf('\n')
 fprintf('The actual sectional lift coefficient is %f (N = 500) \n',CL_a)
 fprintf('\n')
 
 %% Error Analysis
 CL_err = zeros(1,491);
 
for i = 1:491
    CL_err(i) = (CL_a-CL(i))/CL_a * 100;
end
i1 = find(CL_err<=1);
N_a = i1(1);
fprintf('The # of discrete vorticies required for less than one percent error is %0.1f \n',N_a);

% NACA-0006
alphar = linspace(0,pi/4); 
cl = zeros(1,length(alphar));
 t1 = 6/100;
 
        XB_half1 = linspace(c,0,N_a/2);
        XB_half2 = flip(XB_half1(1:length(XB_half1)-1));
        
        YB_half1 = (t1/0.2)*c .* (  0.2969*sqrt(XB_half1/c) - 0.126*(XB_half1/c) - 0.3516*(XB_half1/c).^2 + ...
            0.2843*(XB_half1/c).^3 - 0.1036*(XB_half1/c).^4); 
        YB_half2 = (t1/0.2)*c .* (  0.2969*sqrt(XB_half2/c) - 0.126*(XB_half2/c) - 0.3516*(XB_half2/c).^2 + ...
            0.2843*(XB_half2/c).^3 - 0.1036*(XB_half2/c).^4);
       XB = [XB_half1 XB_half2];
       YB = [-YB_half1 YB_half2];
       
       for i = 1:length(alphar)
    [cl1(i),~,~,~,~] = vortex_panel(XB,YB,V_inf,rad2deg(alphar(i)),0);
       end

 % NACA-0012      
    t2 = 12/100;
     XB_half1 = linspace(c,0,N_a/2);
        XB_half2 = flip(XB_half1(1:length(XB_half1)-1));
        
        YB_half1 = (t2/0.2)*c .* (  0.2969*sqrt(XB_half1/c) - 0.126*(XB_half1/c) - 0.3516*(XB_half1/c).^2 + ...
            0.2843*(XB_half1/c).^3 - 0.1036*(XB_half1/c).^4); 
        YB_half2 = (t2/0.2)*c .* (  0.2969*sqrt(XB_half2/c) - 0.126*(XB_half2/c) - 0.3516*(XB_half2/c).^2 + ...
            0.2843*(XB_half2/c).^3 - 0.1036*(XB_half2/c).^4);
       XB = [XB_half1 XB_half2];
       YB = [-YB_half1 YB_half2];
       
  for i = 1:length(alphar)
    [cl2(i),~,~,~,~] = vortex_panel(XB,YB,V_inf,rad2deg(alphar(i)),0);
end
       
 % NACA-0024   
    t3 = 24/100;
  XB_half1 = linspace(c,0,N_a/2);
        XB_half2 = flip(XB_half1(1:length(XB_half1)-1));
        
        YB_half1 = (t3/0.2)*c .* (  0.2969*sqrt(XB_half1/c) - 0.126*(XB_half1/c) - 0.3516*(XB_half1/c).^2 + ...
            0.2843*(XB_half1/c).^3 - 0.1036*(XB_half1/c).^4); 
        YB_half2 = (t3/0.2)*c .* (  0.2969*sqrt(XB_half2/c) - 0.126*(XB_half2/c) - 0.3516*(XB_half2/c).^2 + ...
            0.2843*(XB_half2/c).^3 - 0.1036*(XB_half2/c).^4);
       XB = [XB_half1 XB_half2];
       YB = [-YB_half1 YB_half2];
 
 
for i = 1:length(alphar)
    [cl3(i),~,~,~,~] = vortex_panel(XB,YB,V_inf,rad2deg(alphar(i)),0);
end

% Plot Together
figure(25)
plot(alphar,cl1)
hold on
plot(alphar,cl2)
hold on
plot(alphar,cl3)
xlabel('Angle of Attack (rad)')
ylabel('Sectional Lift Coefficient')
title('Alpha vs Cl for three different NACA airfoils')
legend('NACA 0006','NACA 0012','NACA0024')
