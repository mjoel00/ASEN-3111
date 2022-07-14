function [e,c_L,c_Di] = pllt(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)

    Nlin = 1:N;
    theta = (Nlin*pi)./(2*N);
% convert to rad 
aero_t = deg2rad(aero_t); 
aero_r = deg2rad(aero_r);
geo_t = deg2rad(geo_t);
geo_r = deg2rad(geo_r);

% Aspect Ratio
AR = (2*b)/(c_r + c_t);

% Get linear eqns 
c = c_r + (c_t - c_r) * cos(theta);
%c = flip(c);
a0 = a0_r + (a0_t - a0_r) * cos(theta);
%a0 = flip(a0);
aero = aero_r + (aero_t - aero_r) * cos(theta);
%aero = flip(aero);
geo = geo_r + (geo_t - geo_r) * cos(theta);
%geo = flip(geo);

% Solve for A_n coefficients

num = zeros(1,N);
den = zeros(N,N);
for i = 1:N
    num(i)=geo(i)-aero(i);
    for j = 1:N
        n = (2*j-1);
        den(i,j) = ((4*b*sin(n*theta(i)))/(a0(i)*c(i))) + ((n*sin(n*theta(i)))/(sin(theta(i)))); 
    end  
end
A =  inv(den) * num';

%%  Outputs
    
    % coefficient of lift
    c_L = A(1)*pi*AR;
    
    % drag factor
    delta = 0;
    for j = 2:N
        n = (2*j-1);
        delta = delta + n*((A(j)/A(1))^2);
    end 
    
    % solve for span efficiency
    e = 1 / (1+delta);
    
    % Solve for induced drag coefficient
    c_Di = c_L^2/(pi*e*AR);
    
end