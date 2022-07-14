%% ASEN 3111 - Computational Assignment 1 - Main
% To find the sectional lift and drag coefficients of a rotating cylinder.
% Further, Produce plots of the sectioonal lift and drag coefficients vs
% the number of panels used to discretize the surface of the cylinder.
% Using both the trapezoidal rule and Simpson's rule. Lastly, find the
% number of panels required to acheive a predicted sectioanl lift
% coefficient with a 1/10 percent relative error using the trapezoidal rule
%
% Author: Cole MacPherson
% Collaborators: R. Block, Z. Lesan, S. Mansfield, A. Uprety
% Date: 26th Jan 2021

%% Housekeeping

tic
clc;
clear;
close all;

%% Finding Sectional Lift and Drag Coefficients

fprintf("PROBLEM 1:\n\nFinding sectional coefficients of lift and drag...\n\n");

syms theta R V_inf 

% housekeeping w/ variables
Gamma = -2*pi*R*V_inf; 
C_p = 1 - ( 4*sin(theta)^2 + ((2*Gamma*sin(theta))/(pi*R*V_inf)) + (Gamma/(2*pi*R*V_inf))^2 );

% calculating sectional coefficients
c_l = (-1/2)*int(C_p*sin(theta),theta,0,2*pi);
c_d = (-1/2)*int(C_p*cos(theta),theta,0,2*pi);

fprintf("The sectional coefficient of lift was fould to be: " + char(c_l) + "\n");
fprintf("The sectional coefficient of drag was fould to be: " + char(c_d) + "\n\n");

%% Sectional Coefficients vs Number of Panels (Trapezoidal Rule)

fprintf("Implementing trapezoidal rule...\n\n");

% initializing variables
a = 0;
b = 2*pi;
N = 1:50;
c_l_trapz = zeros(1,length(N));
c_d_trapz = zeros(1,length(N));

for i = 1:length(N)
    
    % find x_value vector for given panel
    x_values = zeros(1,i+1);
    width = (b-a)/(i);
    for j = 1:(i+1)
       x_values(j) = a + width*(j-1); 
    end
       
    for k = 1:i
        
        % coeff of pressure
        C_p_1_trapz = 4*sin(x_values(k))-4*sin(x_values(k))^2;
        C_p_2_trapz = 4*sin(x_values(k+1))-4*sin(x_values(k+1))^2;
        
        % f(theta)
        f_1_cl_trapz = C_p_1_trapz*sin(x_values(k));
        f_2_cl_trapz = C_p_2_trapz*sin(x_values(k+1));
        f_1_cd_trapz = C_p_1_trapz*cos(x_values(k));
        f_2_cd_trapz = C_p_2_trapz*cos(x_values(k+1));
        
        % coeff of lift
        temp = (-1/2)*(x_values(k+1)-x_values(k))*((f_2_cl_trapz+f_1_cl_trapz)/2);
        c_l_trapz(i) = c_l_trapz(i) + temp;
        
        % coeff of drag
        temp = (-1/2)*(x_values(k+1)-x_values(k))*((f_2_cd_trapz+f_1_cd_trapz)/2);
        c_d_trapz(i) = c_d_trapz(i) + temp;
        
    end
    
end

fprintf("Ploting findings...\n\n");

% plot results
fignum = 1;
figure(fignum)
subplot(2,1,1)
plot(N,c_l_trapz,'o','linewidth',2)
title('Sectional Coefficients vs Number of Intergation Panels (Trapezodial Rule)')
ylabel('c_{l}');
subplot(2,1,2)
plot(N,c_d_trapz,'o','linewidth',2)
ylabel('c_{d}');
xlabel('Number of Intergration Panels');

fignum = fignum + 1;

%% Sectional Coefficients vs Number of Panels (Simpson's Rule)

fprintf("Implementing simpson's rule...\n\n");

% constant declaration
N = 1:50;
c_l_simp = zeros(1,length(N));
c_d_simp = zeros(1,length(N));

for i = 1:length(N)
    
    % find x_value vector for given panel
    x_values = zeros(1,i+1);
    width = (b-a)/(i);
    for j = 0:i
       x_values(j+1) = a + width*(j); 
    end
    h = (b-a)/(2*i);
    
    for k = 1:i
        
        % define sub-panel point
        midpoint = (x_values(k)+x_values(k+1))/2;
        
        % coeff of pressure
        c_p_1_simp = 4*sin(x_values(k))-4*sin(x_values(k))^2;
        c_p_2_simp = 4*sin(midpoint)-4*sin(midpoint)^2;
        c_p_3_simp = 4*sin(x_values(k+1))-4*sin(x_values(k+1))^2;
        
        % f(theta)
        f_1_cl_simp = c_p_1_simp*sin(x_values(k));
        f_2_cl_simp = c_p_2_simp*sin(midpoint);
        f_3_cl_simp = c_p_3_simp*sin(x_values(k+1));
        f_1_cd_simp = c_p_1_simp*cos(x_values(k));
        f_2_cd_simp = c_p_2_simp*cos(midpoint);
        f_3_cd_simp = c_p_3_simp*cos(x_values(k+1));
        
        % coeff of lift
        temp = f_1_cl_simp + 4*f_2_cl_simp + f_3_cl_simp;
        c_l_simp(i) = c_l_simp(i) + temp;
        
        % coeff of drag
        temp = f_1_cd_simp + 4*f_2_cd_simp + f_3_cd_simp;
        c_d_simp(i) = c_d_simp(i) + temp;
        
    end
    
    % multiply by h/3
    c_l_simp(i) = (-1/2)*(c_l_simp(i)*(h/3));
    c_d_simp(i) = (-1/2)*(c_d_simp(i)*(h/3));    
    
end

fprintf("Ploting findings...\n\n");

% plot results
figure(fignum)
subplot(2,1,1)
plot(N,c_l_simp,'o','linewidth',2)
title('Sectional Coefficients vs Number of Intergation Panels (Simpson Rule)')
ylabel('c_{l}');
subplot(2,1,2)
plot(N,c_d_simp,'o','linewidth',2)
ylabel('c_{d}');
xlabel('Number of Intergration Panels');

fignum = fignum + 1;

%% Number of Panels of 1/10 Percent Error

fprintf("Finding the number of panels required for both methods\nto have less than a 1/10 percent error...\n\n");

% initialization
err_trapz = 1;
err_simp = 1;
err = (1/10)/100;
i = 1;
c_l = double(c_l);
c_d = double(c_d);
trapz_bool = false;
simp_bool = false;

while err_trapz > err || err_simp > err
    
    trapz_cl = abs((c_l-c_l_trapz(i))/c_l);
    simp_cl = abs((c_l-c_l_simp(i))/c_l);
    
    if trapz_cl <= err && trapz_bool == false
        err_trapz = trapz_cl;
        trapz_panel = i;
        trapz_bool = true;
    else
        err_trapz = trapz_cl;
    end
    
    if simp_cl <= err && simp_bool == false
        err_trapz = simp_cl;
        simp_panel = i;
        simp_bool = true;
    else
        err_simp = simp_cl;
    end
    
    i = i + 1;
    
end

fprintf("The number of panels required for the trapezoidal rule\n to have less than 1/10 percent error is %i \n",trapz_panel);
fprintf("The number of panels required for the simpson's rule\n to have less than 1/10 percent error is %i \n\n",simp_panel);

%% Find Lift and Drag From Cp Data

fprintf("PROBLEM 2:\n\nFinding the lift and drag using Cp data...\n\n");

load("Cp.mat"); % load given file

% intialization variables
c = 5; % chord length [m]
AoA = 9; % angle of attack [degrees]
V_inf = 20; % freestream velocity [m/s]
rho_inf = 1.225; % freestream density [kg/m^3]
P_inf = 101.3e3; % freestream pressure [Pa]
q_inf = (1/2)*rho_inf*V_inf^2; % dynamic pressure [Pa]
t = 12/100; % thickness to chord ratio

% initialization of key vectors used for analysis
N = 1:1000;
N_prime_u = zeros(1,length(N));
N_prime_l = zeros(1,length(N));
N_prime = zeros(1,length(N));
A_prime_u = zeros(1,length(N));
A_prime_l = zeros(1,length(N));
A_prime = zeros(1,length(N));
L_prime = zeros(1,length(N));
D_prime = zeros(1,length(N));

for i = 1:length(N)
    
    % find x_value vector for given panel laong with teh panel width
    x_values = zeros(1,i+1);
    width = c/length(N);
    for j = 0:i
       x_values(j+1) = width*(j); 
    end
    
    % get pressure for the upper and lower portions of the airfoil
    P_u = (fnval(Cp_upper, x_values/c)*q_inf)+P_inf;
    P_l = (fnval(Cp_lower, x_values/c)*q_inf)+P_inf;
    
    for k = 1:i-1
        
        % find the sectional normal force of upper and lower surfaces of
        % the airfoil
        N_prime_u(i) = N_prime_u(i) + (width/2)*(P_u(k+1)+P_u(k));
        N_prime_l(i) = N_prime_l(i) + (width/2)*(P_l(k+1)+P_l(k));
        
        % find the y values of the airfoil at point x, to solve for the
        % sectional axial force
        y_t_1 = (t/0.2)*c*(0.2969*sqrt(x_values(k)/c)-0.126*(x_values(k)/c)...
            -0.3516*(x_values(k)/c)^2+0.2843*(x_values(k)/c)^3-0.1036*(x_values(k)/c)^4);
        y_t_2 = (t/0.2)*c*(0.2969*sqrt(x_values(k+1)/c)-0.126*(x_values(k+1)/c)...
            -0.3516*(x_values(k+1)/c)^2+0.2843*(x_values(k+1)/c)^3-0.1036*(x_values(k+1)/c)^4);
        
        % find the sectional axial force of upper and lower surfaces of
        % the airfoil
        A_prime_u(i) = A_prime_u(i) + ((y_t_2-y_t_1)/2)*(-P_u(k+1)-P_u(k));
        A_prime_l(i) = A_prime_l(i) + ((y_t_2-y_t_1)/2)*(P_l(k+1)+P_l(k));
        
    end
    
    % consolidate the upper and lower surface vectors into one vector for
    % the sectional normal and axial forces
    N_prime(i) = N_prime_l(i) - N_prime_u(i);
    A_prime(i) = A_prime_l(i) + A_prime_u(i);
    
    % solving for lift and drag with consolidated sectional normal and
    % axial force vectors
    L_prime(i) = N_prime(i)*cosd(AoA) - A_prime(i)*sind(AoA);
    D_prime(i) = N_prime(i)*sind(AoA) - A_prime(i)*cosd(AoA);
    
end

fprintf("The sectional lift was found to be %0.2f [N/m]\n",L_prime(end));
fprintf("The sectional drag was found to be %0.2f [N/m]\n\n",D_prime(end));

%% Error Integration Points

fprintf("Find the amount of integration points rquired for a 5%%, 1%%, and 1/10%% relative error...\n\n");

% initialization of variables for error analysis
percent5 = false;
point5 = 0;
percent1 = false;
point1 = 0;
percent1tenth = false;
point1thenth = 0;
i = 1;

while percent5 == false || percent1 == false || percent1tenth == false
    
    L_err = abs((L_prime(end)-L_prime(i))/L_prime(end)); % error between point and actual
        
    if L_err <= 0.05 && percent5 == false
        fprintf("The number of integration points required for a 5%% error is %i points\n",i);
        percent5 = true;
    end
    
    if L_err <= 0.01 && percent1 == false
        fprintf("The number of integration points required for a 1%% error is %i points\n",i);
        percent1 = true;
    end
    
    if L_err <= (1/10)/100 && percent1tenth == false
        fprintf("The number of integration points required for a 1/10%% error is %i points\n\n",i);
        percent1tenth = true;
    end
    
    i = i + 1;
    
end