% Matthew J. Pabin
% ASEN 3111 - Aerodynamics 
% CA-01 : Computation of Lift and Drag
% Submitted Jan 30 2022

clear
clc
close all

%% Problem 1
syms theta

%1A
c_l = @(theta) sin(theta).*(1 - (4*sin(theta).^2 + 4*sin(theta) + 1));
c_d = @(theta) cos(theta).*(1 - (4*sin(theta).^2 + 4*sin(theta) + 1));

c_l = 0.5*integral(c_l,0,2*pi); 
c_d = 0.5*integral(c_d,0,2*pi);

fprintf('-------Question 1--------')
fprintf('\n')
fprintf('C_L = %0.4f = -2\x03C0\n',c_l)
fprintf('\n')
fprintf('C_D = %0.0f \n',c_d)
fprintf('\n')


%1B

N = 1000;
r = linspace(0,2*pi,N);
c_l_T = zeros(1,N);
c_d_T = zeros(1,N);

% Coefficient of Lift Calculation using trap rule
for  i=1:N-1
    c_l_T(i) = 0.5*(r(i+1) - r(i)).* ...
    (sin(r(i+1)).*(1 - (4*sin(r(i+1)).^2 + ...
    4*sin(r(i+1)) + 1)) + sin(r(i)).* ...
    (1 - (4*sin(r(i)).^2 + 4*sin(r(i)) + 1)))/2;
end

% Coefficient of Drag Calculation using trap rule
for  i=1:N-1
    c_d_T(i) = 0.5*(r(i+1) - r(i)).* ...
    (cos(r(i+1)).*(1 - (4*sin(r(i+1)).^2 + ...
    4*sin(r(i+1)) + 1)) + cos(r(i)).* ...
    (1 - (4*sin(r(i)).^2 + 4*sin(r(i)) + 1)))/2;
end

for i = 1:(length(c_l_T))
   clTSUM(i) = sum(c_l_T(1:i));
end

for i = 1:(length(c_d_T))
   cdTSUM(i) = sum(c_d_T(1:i));
end

figure(1)
plot(1:N,clTSUM,'b','LineWidth',2)
title('Trapezoidal Rule')
xlabel('# of Iterations')
ylabel('Coefficient of Lift')
ylim([-8 1])
figure(2)
plot(1:N,cdTSUM,'r','LineWidth',2)
title('Trapezoidal Rule')
xlabel('# of Iterations')
ylabel('Coefficient of Drag')
ylim([-8 1])

% 1C
c_l_S = zeros(1,N);
c_d_S = zeros(1,N);

% Coef of drag using simpsons
for  i=2:N-1
    c_l_S(i) = ((2*pi) / (6*N)).*(sin(r(i-1)).* ...
    (1 - (4*sin(r(i-1)).^2 + 4*sin(r(i-1)) + 1)) + ...
    4*sin(r(i)).* ...
    (1 - (4*sin(r(i)).^2 + 4*sin(r(i)) + 1)) + ...
    sin(r(i+1)).*(1 - (4*sin(r(i+1)).^2 + 4*sin(r(i+1)) + 1)))/2;
end

% Coef of Drag using Simpsons 
for  i=2:N-1
    c_d_S(i) = ((2*pi) / (6*N)).*(cos(r(i-1)).* ...
    (1 - (4*sin(r(i-1)).^2 + 4*sin(r(i-1)) + 1)) + ...
    4*cos(r(i)).* ...
    (1 - (4*sin(r(i)).^2 + 4*sin(r(i)) + 1)) + ...
    cos(r(i+1)).*(1 - (4*sin(r(i+1)).^2 + 4*sin(r(i+1)) + 1)))/2;
end


for i = 1:(length(c_l_S))
   clSSUM(i) = sum(c_l_S(1:i));
end

for i = 1:(length(c_d_S))
   cdSSUM(i) = sum(c_d_S(1:i));
end


figure(3)
plot(1:N,clSSUM,'b','LineWidth',2)
title('Simpsons Rule')
xlabel('# of Iterations')
ylabel('Coefficient of Lift')
ylim([-8 1])
figure(4)
plot(1:N,cdSSUM,'r','LineWidth',2)
title('Simpsons Rule')
xlabel('# of Iterations')
ylabel('Coefficient of Drag')
ylim([-8 1])


% 1D
c_l_err = c_l * 0.01;
for i = 1:N
if clTSUM(i) <= (c_l - c_l_err) && clTSUM(i) >= (c_l + c_l_err)
    fprintf('# of integration Points Needed for less than 1 percent relative error using Trapezoidal Rule : %0.0f \n',i)
    break
end
end
fprintf('\n')
for i = 1:N
if clSSUM(i) <= (c_l - c_l_err) && clSSUM(i) >= (c_l + c_l_err)
    fprintf('# of integration Points Needed for less than 1 percent relative error using Simpson''s Rule : %0.0f \n',i)
    break
end
end
fprintf('\n')
%% Problem 2

% import Cp.mat
load('Cp.mat');

N = 1000;
c = 2; % 2m
alpha = deg2rad(9); % rad
V_inf = 60; % m/s
rho_inf = 1; % kg/m^3
p_inf = 85.5*10^3; % Pa
t = 12/100; % Chord thickness

% Upper Surface
yu = @(x) (t*c)/0.2 * (0.2969 * sqrt(x./c) - 0.1260 * (x./c) ...
   - 0.3516 * (x./c).^2 + 0.2843 * (x./c).^3 - 0.1036 * (x./c).^4);

r = linspace(0,c,N);
delta_yu = double(yu(r));

Cpu = fnval(Cp_upper,r/c);

% Lower Surface
yl = @(x) (-t*c)/0.2 * (0.2969 * sqrt(x./c) - 0.1260 * (x./c) ...
    - 0.3516 * (x./c).^2 + 0.2843 * (x./c).^3 - 0.1036 * (x./c).^4);

delta_yl = double(yl(r));

Cpl = fnval(Cp_lower,r/c);

% loop thru iterations
for k = 1:1000
% Normal/Axial Forces
for i = 1:k-1
    
   Nl(i) = (r(i+1) - r(i)) * (Cpl(i+1) ...
                + Cpl(i))/2;
   Nu(i) = (r(i+1) - r(i)) * (Cpu(i+1) ...
                + Cpu(i))/2;
    
   Al(i) = ((r(i+1) - r(i)) * (Cpl(i+1) ...
                + Cpl(i))/2) * (-(delta_yl(i) ...
                - delta_yl(i+1)) / (r(i+1) - r(i)));
   Au(i) = ((r(i+1) - r(i)) * (Cpu(i+1) ...
                + Cpu(i))/2) * (-(delta_yu(i) ...
                - delta_yu(i+1)) / (r(i+1) - r(i)));
    
end


if k == 1
    Cn(k) = 0;
    Ca(k) = 0;
else
    Cn(k) = (1/c) * (sum(Nl) - sum(Nu));
    Ca(k) = (1/c) * (sum(Au) - sum(Al));
end

% Coef of lift and drag
Cl = Cn * cos(alpha) - Ca * sin(alpha);
Cd = Cn * sin(alpha) + Ca * cos(alpha);
% Lift
Lift = 0.5 * rho_inf * V_inf^2 * c * Cl;
Drag = 0.5 * rho_inf * V_inf^2 * c * Cd;
end

% Base values for lift and drag
L = Lift(end);
D = Drag(end);
fprintf('\n')
fprintf('\n')
fprintf('-------Question 2--------\n')
fprintf('L = %0.4f (N/m)\n',L)
fprintf('\n')
fprintf('D = %0.4f (N/m)\n',D)
fprintf('\n')

% Error 
Cl_err = abs(Cl - Cl(end)) ./ Cl(end);

% Vectors of values within acceptable tolerance
x1 = find(Cl_err < 0.1);
x2 = find(Cl_err < 0.01);
x3 = find(Cl_err < 0.001);

% Print first value of each vector 
fprintf('# of integration points needed for less than 10 percent relative error:  %0.0f \n',x1(1))
fprintf('\n')
fprintf('# of integration points needed for less than 1 percent relative error:  %0.0f \n',x2(1))
fprintf('\n')
fprintf('# of integration points needed for less than 0.1 percent relative error:  %0.0f \n',x3(1))
fprintf('\n')
 
% end



