%% ASEN 3111
% Computational Assignment 3: Flow over Finite Wings
% Matthew J. Pabin
% 3/7/22

clear
clc
close all

% Givens

V =  60;   %m/s
c =  2;    %m
N =  139;     % From CA2

% NACA 0012
m00 = 0/100;
p00 = 0/10;
t00 = 12/100;

% NACA 2412
m24 = 2/100;
p24 = 4/10;
t24 = 12/100;

% NACA 4412
m44 = 4/100;
p44 = 4/10;
t44 = 12/100;

 % Problem 1

% call NACA_airfoils
[x00,y00,aLo00,dCLda00] = NACA_Airfoils(m00,p00,t00,c,N);
[x24,y24,aLo24,dCLda24] = NACA_Airfoils(m24,p24,t24,c,N);
[x44,y44,aLo44,dCLda44] = NACA_Airfoils(m44,p44,t44,c,N);
%
alpha = linspace(-10,75,N);
Cl00t = dCLda00 * (deg2rad(alpha) - aLo00);
Cl24t = dCLda24 * (deg2rad(alpha) - aLo24);
Cl44t = dCLda44 * (deg2rad(alpha) - aLo44);
Cl00 = zeros(1,length(alpha));
Cl24 = zeros(1,length(alpha));
Cl44 = zeros(1,length(alpha));

for i = 1:length(alpha)
[Cl00(i),~,~,~,~] = vortex_panel(x00,y00,V,alpha(i),0);
[Cl24(i),~,~,~,~] = vortex_panel(x24,y24,V,alpha(i),0);
[Cl44(i),~,~,~,~] = vortex_panel(x44,y44,V,alpha(i),0);
end

X = [ dCLda00 aLo00; dCLda24 aLo24; dCLda44 aLo44];
X = array2table(X,'VariableNames',{'Lift Slope (1/rad)','Zero-Lift Angle of Attack (rad)'},'RowNames',...
    {'NACA 0012','NACA 2412','NACA 4412'});

fprintf('----- Problem 1 -----')
fprintf('\n')
fprintf('\n')
disp(X);

figure(1);
plot(alpha,Cl00,'r')
hold on
plot(alpha,Cl24,'b')
hold on
plot(alpha,Cl44,'g')
hold on 
plot(alpha,Cl00t,'r--')
hold on
plot(alpha,Cl24t,'b--')
hold on
plot(alpha,Cl44t,'g--')
title('\alpha vs. C_{L} ')
legend('NACA 0012','NACA 2412','NACA 4412')
xlabel('\alpha (deg)')
ylabel('C_{L}')



%% Problem 2

% Givens
N = 100;
V = 80 * 1.688; % ft/s 
rho = 0.001756; % slug/ft^3
b = 33.333; % ft
c_r = 5.333; %ft
c_t = 3.708333; %ft
a0_r = dCLda24;  % 1/rad
a0_t = dCLda00;  % 1/rad
aero_r = rad2deg(aLo24); %deg
aero_t = rad2deg(aLo00); %deg
geo_r = 1;  %deg
geo_t = 0;  %deg

% call PLLT
[e,c_L,c_Di] = pllt(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

S = 2*((0.5*b*c_t) + (0.25*b*(c_r-c_t)));
q = 0.5*rho*V^2;
L = c_L * q * S;
Di = c_Di * q * S;

c_Lvec = zeros(1,N);
for l = 1:N
    [~,c_Lvec(l),~] = pllt(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,l);
end 
err1 = c_L*.1;
err2 = c_L*.01;
err3 = c_L*0.001;

c_Lerr = (c_Lvec - c_L);
N1 = find(c_Lerr < err1);
N1 = N1(1);
N2 = find(c_Lerr < err2);
N2 = N2(1);
N3 = find(c_Lerr < err3);
N3 = N3(1);

fprintf('----- Problem 2 -----')
fprintf('\n')
fprintf('\n')
fprintf('Lift: <strong> %f </strong> (lbf)',L)
fprintf('\n')
fprintf('\n')
fprintf('Induced Drag: <strong> %f </strong> (lbf)',Di)
fprintf('\n')
fprintf('\n')
fprintf('# of Odd Terms required for solution with < 10 percent error: <strong> %0.f </strong>',N1);
fprintf('\n')
fprintf('# of Odd Terms required for solution with < 1 percent error: <strong> %0.f </strong>',N2);
fprintf('\n')
fprintf('# of Odd Terms required for solution with < 0.1 percent error: <strong> %0.f </strong>',N3);
fprintf('\n')
fprintf('\n')

%% Problem 3

fprintf('----- Problem 3 -----')
fprintf('\n')
fprintf('\n')

% call pllt again to find e 
N = 70;
% taper ratio vector
taprat = linspace(0,1,N);
% initialize AR
AR = [4 6 8 10];
c_r = 1;
% Loop through each taper ratio for each AR
for i = 1:length(AR)
    for j = 1:N
        c_t = c_r * taprat(j);
        b = 0.5*(c_t + c_r)*AR(i);
        [e(i,j),~,~] = pllt(b,a0_t,a0_r,c_t,c_r,0,0,1,1,N);
    end
    % plot
    figure(2);
    plot(taprat,e(i,:))
    hold on
end
    title('Span Efficiency Factor vs Taper Ratio')
    xlabel('c_{t}/c_{r}')
    ylabel('e')
    legend('AR = 4','AR = 6','AR = 8','AR = 10')
    grid on;
    
    
