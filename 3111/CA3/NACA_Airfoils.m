function [x,y,aLo,dCLda] = NACA_Airfoils(m,p,t,c,N)
% Boundary Points without mean camber line
        
     xb = linspace(0,c,N);
        
        
     yb = (t/0.2)*c .* (  0.2969*sqrt(xb/c) - 0.126*(xb/c) - 0.3516*(xb/c).^2 + ...
            0.2843*(xb/c).^3 - 0.1036*(xb/c).^4); 
        
     
% Mean Camber line 
z = zeros(1,length(xb));
if (m==0) && (p==0)
z = zeros(1,length(xb));
else
    for i = 1:length(xb)
        if xb(i) <=  p*c
            z(i) = m .* (xb(i)./p^2) .* (2*p - (xb(i)./c));
        elseif xb(i) > p*c    
            z(i) = m .* ((c-xb(i))./(1-p)^2) .* (1 + (xb(i)./c) - 2*p);     
        end
     end       
end
 
% Derivative of mean camber line
dzdx = zeros(1,length(xb));
if (m==0) && (p==0)
    dzdx = zeros(1,length(xb));
else
    for i = 1:length(xb)
        if xb(i) <= p*c
        dzdx(i) = 2*m.*(c*p-xb(i)) ./ (c*p^2);
        elseif xb(i) > p*c
        dzdx(i) = 2*m.*(c*p-xb(i)) ./ (c*(1-p)^2);
        end
    end
end
% Solve for zeta
zeta = atan(dzdx);
% Solve for airfoil coordinates w mean camber line
xu = xb - yb.*sin(zeta);
xl = xb + yb.*sin(zeta);
yu = z + yb.*cos(zeta);
yl = z - yb.*cos(zeta);

x = [flip(xl),xu(2:end)];
y = [flip(yl),yu(2:end)];

%Change of vars
theta = linspace(0,pi,N);
dzdtheta = zeros(1,length(xb));
if (m==0) && (p==0)
    dzdtheta = zeros(1,length(xb));
else
    for i = 1:length(xb)
        if xb(i) <= p*pi
    dzdtheta(i) = (m * (-1 + 2*p + cos(theta(i))))/p^2;
        elseif xb(i) >= p*pi
    dzdtheta(i) = (m * (-1 + 2*p + cos(theta(i))))/(-1 + p)^2;
        end
    end    
end

aLo = -(1/pi) .* trapz(theta,dzdtheta.*(cos(theta)-1));
dCLda = 2*pi;

end