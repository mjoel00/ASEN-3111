% The purpose of this computational lab assignment is to complete a basic
% numerical integration example on a given airfoil to form a foundation
% going forward. In this lab we used Trapezoidal Rule Method and Simpson's
% Rule Method to numerically integrate over a stationary NACA 0012 airfoil
% to understand how both methods calculate the coefficients of drag and
% lift, as well as the actual lift and drag forces.
%
%   Author: Ajay Dhindsa
%   Collaborators: Abdullah Almugairin
%   Date: 17 September 2020

clear all
clc
close all

%% Question 1A - Analytic Solution of c_l and c_d

% Creating C_l and C_d functions to be integrated
syms theta

c_l_a = sin(theta).*(1 - (4*sin(theta).^2 + 4*sin(theta) + 1));
c_d_a = cos(theta).*(1 - (4*sin(theta).^2 + 4*sin(theta) + 1));

% Integrating C_l and C_d functions analytically
c_l_A = int(c_l_a,theta);
c_d_A = int(c_d_a,theta);

c_l = @(theta) sin(theta).*(1 - (4*sin(theta).^2 + 4*sin(theta) + 1));
c_d = @(theta) cos(theta).*(1 - (4*sin(theta).^2 + 4*sin(theta) + 1));

% Analytic Solutions
c_l_A_E = -0.5*integral(c_l,0,2*pi); % Equals 2pi
c_d_A_E = -0.5*integral(c_d,0,2*pi); % Equals 0

% Displaying Results
fprintf('Question 1A Results: ')
fprintf('\n')
fprintf('\n')

fprintf('Analytic Solution for c_l = 2\x03C0 %g\t')
fprintf('\n')
fprintf('Analytic Solution for c_d = 0')
fprintf('\n')
fprintf('\n')
fprintf('\n')

%% Question 1B - Composite Trapezoidal Rule

% Setting number of panels and base for trap integration method
N = 1000;
range = linspace(0,2*pi,N);
c_l_T = zeros(1,N);

% Coefficient of Lift Calculation using trap rule
for  i=1:N-1
    c_l_T(i) = 0.5*(range(i+1) - range(i)).* ...
    (sin(range(i+1)).*(1 - (4*sin(range(i+1)).^2 + ...
    4*sin(range(i+1)) + 1)) + sin(range(i)).* ...
    (1 - (4*sin(range(i)).^2 + 4*sin(range(i)) + 1)));
end

% Coefficient of Drag Calculation using trap rule
for  i=1:N-1
    c_d_T(i) = 0.5*(range(i+1) - range(i)).* ...
    (cos(range(i+1)).*(1 - (4*sin(range(i+1)).^2 + ...
    4*sin(range(i+1)) + 1)) + cos(range(i)).* ...
    (1 - (4*sin(range(i)).^2 + 4*sin(range(i)) + 1)));
end

% Multiplying final -1/2 in front of integral
c_l_T = c_l_T*-0.5;
c_d_T = c_d_T*-0.5;

% Creating sum vectors to plot against number of panels
for i = 1:(length(c_l_T)-1)
   c_l_T_SUM(i) = sum(c_l_T(1:i));
end

for i = 1:(length(c_d_T)-1)
   c_d_T_SUM(i) = sum(c_d_T(1:i));
end

% Plotting Results
figure(1)
plot(1:length(c_l_T_SUM),c_l_T_SUM)
title('Coefficient of Lift, $C_L$ VS Number of Panels, $N$ - Composite Trapezoidal Rule', 'Interpreter', 'latex')
xlabel('Number of Panels, $N$','Interpreter','latex')
ylabel('Coefficient of Lift, $C_L$', 'Interpreter', 'latex')

figure(2)
plot(1:length(c_d_T_SUM),c_d_T_SUM)
title('Coefficient of Drag, $C_D$ VS Number of Panels, $N$ - Composite Trapezoidal Rule', 'Interpreter', 'latex')
xlabel('Number of Panels, $N$','Interpreter','latex')
ylabel('Coefficient of Drag, $C_D$', 'Interpreter', 'latex')

%% Question 1C - Composite Simpson's Rule

% Number of panels and base for simpsons rule
N = 1000;
range = linspace(0,2*pi,N);

% Coefficient of Lift calculation using simpson's rule
for  i=2:N-1
    c_l_S(i) = ((2*pi) / (6*N)).*(sin(range(i-1)).* ...
    (1 - (4*sin(range(i-1)).^2 + 4*sin(range(i-1)) + 1)) + ...
    4*sin(range(i)).* ...
    (1 - (4*sin(range(i)).^2 + 4*sin(range(i)) + 1)) + ...
    sin(range(i+1)).*(1 - (4*sin(range(i+1)).^2 + 4*sin(range(i+1)) + 1)));
end

% Coefficient of Drag calculation using simpson's rule
for  i=2:N-1
    c_d_S(i) = ((2*pi) / (6*N)).*(cos(range(i-1)).* ...
    (1 - (4*sin(range(i-1)).^2 + 4*sin(range(i-1)) + 1)) + ...
    4*cos(range(i)).* ...
    (1 - (4*sin(range(i)).^2 + 4*sin(range(i)) + 1)) + ...
    cos(range(i+1)).*(1 - (4*sin(range(i+1)).^2 + 4*sin(range(i+1)) + 1)));
end

% Multiplying by -1/2 in front of integral
c_l_S = c_l_S*-0.5;
c_d_S = c_d_S*-0.5;

% Creating sum vectors to plot against number of panels
for i = 1:(length(c_l_S)-1)
   c_l_S_SUM(i) = sum(c_l_S(1:i));
end

for i = 1:(length(c_d_S)-1)
   c_d_S_SUM(i) = sum(c_d_S(1:i));
end

% Plotting Results
figure(3)
plot(1:length(c_l_S_SUM),c_l_S_SUM)
title('Coefficient of Lift, $C_L$ VS Number of Panels, $N$ - Composite Simpson''s Rule', 'Interpreter', 'latex')
xlabel('Number of Panels, $N$','Interpreter','latex')
ylabel('Coefficient of Lift, $C_L$', 'Interpreter', 'latex')

figure(4)
plot(1:length(c_d_S_SUM),c_d_S_SUM)
title('Coefficient of Drag, $C_D$ VS Number of Panels, $N$ - Composite Simpson''s Rule', 'Interpreter', 'latex')
xlabel('Number of Panels, $N$','Interpreter','latex')
ylabel('Coefficient of Drag, $C_D$', 'Interpreter', 'latex')

%% Question 1D

% Finding number of panels for error bound of 0.001 percent
c_l_T_SUM_Error = find((abs(c_l_T_SUM - 2*pi)./(2*pi)) <= 0.0011,1);
c_d_T_SUM_Error = find((abs(c_d_T_SUM - 0.00025)./(0.001)) <= 0.001,1);

% Outputting results for 1D to command window
fprintf('Question 1D Results:');
fprintf('\n')
fprintf('\n')
fprintf(['Number of Panels Required to Obtain Predicted Sectional Lift Coefficient with 0.01 Percent Relative Error Using Composite Trapezoidal Rule: ' num2str(c_l_T_SUM_Error)]);
fprintf('\n')
fprintf(['Number of Panels Required to Obtain Predicted Sectional Drag Coefficient with 0.01 Percent Relative Error Using Composite Trapezoidal Rule: ' num2str(c_d_T_SUM_Error)]);
fprintf('\n')
fprintf('\n')
fprintf('\n')

%% Question 1E

% Finding number of panels for error bound of 0.001 percent
c_l_S_SUM_Error = find((abs(c_l_S_SUM - 2*pi)./(2*pi)) <= 0.0011,1);
c_d_S_SUM_Error = find((abs(c_d_S_SUM - 0.001)./(0.01)) <= 0.001,1);

% Outputting results for 1E to command window
fprintf('Question 1E Results:');
fprintf('\n')
fprintf('\n')
fprintf(['Number of Panels Required to Obtain Predicted Sectional Lift Coefficient with 0.01 Percent Relative Error Using Composite Simpson''s Rule: ' num2str(c_l_S_SUM_Error)]);
fprintf('\n')
fprintf(['Number of Panels Required to Obtain Predicted Sectional Drag Coefficient with 0.01 Percent Relative Error Using Composite Simpson''s Rule: ' num2str(c_d_S_SUM_Error)]);
fprintf('\n')
fprintf('\n')
fprintf('\n')
%% Question 2
% Finding Lift and Drag per unit span of Stationary NACA 0012 Airfoil with
% specified characteristics

% Airfoil Constants
N = 5000; % Number of Panels
c = 2; % 2m
alpha = 9; % aoa of 9 degrees
V_inf = 30; % m/s
rho_inf = 1.225; % kg/m^3
p_inf = 101.3*10^3; % Pa
t = 0.12; % Chord thickness

load Cp

%% Upper Surface
% Upper surface of airfoil calculations
y_upper = @(x) (t*c)/0.2 * (0.2969 * sqrt(x./c) - 0.1260 * (x./c) ...
               - 0.3516 * (x./c).^2 + 0.2843 * (x./c).^3 - 0.1036 * (x./c).^4);

range = linspace(0,c,N);
Delta_y_upper = double(y_upper(range));

Cp_upper_eval = fnval(Cp_upper,range/c);

%% Lower Surface
% Lower surface of airfoil calculations
y_lower = @(x) (-t*c)/0.2 * (0.2969 * sqrt(x./c) - 0.1260 * (x./c) ...
               - 0.3516 * (x./c).^2 + 0.2843 * (x./c).^3 - 0.1036 * (x./c).^4);

Delta_y_lower = double(y_lower(range));

Cp_lower_eval = fnval(Cp_lower,range/c);

%% Trapezoidal Rule
% Doing trapezoidal rule method to integrate both sides of airfoil to
% determine the coefficients of lift and drag

for N = 1:5000

for i = 1:N-1
    
   lower_n(i) = (range(i+1) - range(i)) * (Cp_lower_eval(i+1) ...
                + Cp_lower_eval(i))/2;
   upper_n(i) = (range(i+1) - range(i)) * (Cp_upper_eval(i+1) ...
                + Cp_upper_eval(i))/2;
    
   lower_a(i) = ((range(i+1) - range(i)) * (Cp_lower_eval(i+1) ...
                + Cp_lower_eval(i))/2) * (-(Delta_y_lower(i) ...
                - Delta_y_lower(i+1)) / (range(i+1) - range(i)));
   upper_a(i) = ((range(i+1) - range(i)) * (Cp_upper_eval(i+1) ...
                + Cp_upper_eval(i))/2) * (-(Delta_y_upper(i) ...
                - Delta_y_upper(i+1)) / (range(i+1) - range(i)));
    
end

% Empty value at N = 1 so replaced with 0
if N == 1
    C_n(N) = 0;
    C_a(N) = 0;
else
    C_n(N) = (1/c) * (sum(lower_n) - sum(upper_n));
    C_a(N) = (1/c) * (sum(upper_a) - sum(lower_a));
end

% Final C_l and C_d calculations using C_n and C_a calculated earlier
C_l = C_n * cosd(alpha) - C_a * sind(alpha);
C_d = C_n * sind(alpha) + C_a * cosd(alpha);

Lift = 0.5 * rho_inf * V_inf^2 * c * C_l;

end

Drag = 0.5 * rho_inf * V_inf^2 * c * C_d(end);

% Finding error bounds for C_l and Lift
C_l_Error = abs(C_l - C_l(end)) ./ C_l(end);
Lift_Error = abs(Lift - Lift(end)) ./ Lift(end);

%% Results
% Displaying results in command window
fprintf('Question 2 Results:');
fprintf('\n')
fprintf('\n')

fprintf(['Lift of Stationary NACA 0012: ' num2str(Lift(end)), ' N'])
fprintf('\n')
fprintf(['Drag of Stationary NACA 0012: ' num2str(Drag), ' N'])
fprintf('\n')
fprintf('\n')

fprintf(['Number of Equispaced Integration Points, N, Required to Obtain Lift Solution with 5 Percent Error: ' num2str(find(C_l_Error < 0.05,1))]);
fprintf('\n')

fprintf(['Number of Equispaced Integration Points, N, Required to Obtain Lift Solution with 1 Percent Error: ' num2str(find(C_l_Error < 0.01,1))]);

fprintf('\n')
fprintf(['Number of Equispaced Integration Points, N, Required to Obtain Lift Solution with 0.1 Percent Error: ' num2str(find(C_l_Error < 0.001,1))]);
fprintf('\n')