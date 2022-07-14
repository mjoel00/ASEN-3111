function [c_l,c_dw] = DiamondAirfoil(M, alpha, epsilon1, epsilon2)

gamma = 7/5;



%%%%%% Oblique Shock to Region 1

theta1 = epsilon1 - alpha;

% Wave Angle
beta1 = ObliqueShockBeta(M,theta1,gamma,'Weak');

Mn = M*sind(beta1);

p1_p = 1 + ((2*gamma)/(gamma+1))*(Mn^2 - 1);

Mn1 = sqrt( (1+((gamma-1)/2)*Mn^2) / ((gamma*Mn^2)-((gamma-1)/2)) );

M1 = Mn1 / sind(beta1-theta1);




%%%%%% Oblique Shock to Region 3

theta3 = epsilon1 + alpha;

% Wave Angle
beta3 = ObliqueShockBeta(M,theta3,gamma,'Weak');

Mn = M*sind(beta3);

p3_p = 1 + (((2*gamma)/(gamma+1))*(Mn^2 - 1));

Mn3 = sqrt( (1+((gamma-1)/2)*Mn^2) / ((gamma*Mn^2)-((gamma-1)/2)) );

M3 = Mn3 / sind(beta3-theta3);




%%%%%% Expansion Fan to region 2
theta2 =  (epsilon1 + epsilon2);

[~, v1, ~] = flowprandtlmeyer(gamma,M1, 'mach');

v2 = v1 + theta2;

[M2, ~, ~] = flowprandtlmeyer(gamma,v2, 'nu');

T2_T1 = (1+(((gamma-1)/2)*M1^2)) / (1+(((gamma-1)/2)*M2^2));

p2_p1 = T2_T1^(gamma/(gamma-1));
p2_p = p2_p1 * p1_p;




%%%%%%  Expansion fan to region 4
theta4 =  (epsilon1 + epsilon2);

[~, v3, ~] = flowprandtlmeyer(gamma,M3, 'mach');

v4 = v3 + theta4;

[M4, ~, ~] = flowprandtlmeyer(gamma,v4, 'nu');

T4_T3 = (1+(((gamma-1)/2)*M3^2)) / (1+(((gamma-1)/2)*M4^2));

p4_p3 = T4_T3^(gamma/(gamma-1));
p4_p = p4_p3 * p3_p;




% Calculate Coefficients

% Derive and calculate Cn and Ca due to two different half angles

c_a = (2*sind(epsilon1)*sind(epsilon2))/(gamma*M^2*sind(epsilon1 + epsilon2))* ...
    ((p1_p - p2_p + p3_p - p4_p));

c_n = (-2/(gamma*M^2*sind(epsilon1+epsilon2)))*(((p1_p - p3_p)*cosd(epsilon1)*sind(epsilon2))...
    + ((p2_p - p4_p)*cosd(epsilon2)*sind(epsilon1)));

% calculate Cl and Cdw
c_l = c_n*cosd(alpha) - c_a*sind(alpha);
c_dw = c_n*sind(alpha) + c_a*cosd(alpha);

end