% ASEN 3111: Computational Assignment #4
%
% Matthew J. Pabin
% 4/11/22
% 


% housekeeping
clear
clc 
close all

%% Problem 1

fprintf('----- Problem 1 -----')
fprintf('\n')
fprintf('\n')
fprintf('(Figure 1)')
fprintf('\n')
fprintf('\n')

% Initialize Mach vector

M = [1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.6 1.7 1.8 1.9 2 2.2 2.4 2.6...
    2.8 3.0 3.2 3.4 3.6 3.8 4.0 4.5 5 6 8 10 20 ]; 

% Initialize Deflection angle vector
Thetad = 0:0.1:45.5;
% Ratio of Specific Heats
Gamma = 7/5;
% Is shock wave weak or strong
weak = 'Weak';
strong = 'Strong';

% Initialize beta vectors 
BetadWeak = zeros(length(Thetad),length(M));
BetadStrong = zeros(length(Thetad),length(M));
m = 0;
n = 0;

% Loop through each Mach number and each Theta
for j = 1:length(M)
    for i = 1:length(Thetad)
        
      BetadWeak(i,j) = ObliqueShockBeta(M(j),Thetad(i),Gamma,weak);
      BetadStrong(i,j) = ObliqueShockBeta(M(j),Thetad(i),Gamma,strong);
      
        if  (BetadWeak(i,j) >= 0) && (isreal(BetadWeak(i,j))) 
            m = m + 1;
            tempBetaWeak(m,j) = BetadWeak(i,j);
            tempThetaWeak(m,j) = Thetad(i);
        end
        
        if  (BetadStrong(i,j) >= 0) && (isreal(BetadStrong(i,j))) 
            n = n + 1;
            tempBetaStrong(n,j) = BetadStrong(i,j);
            tempThetaStrong(n,j) = Thetad(i);
        end
        
     
        
    end
    % Plot both weak and strong
    figure(1);
     plot([tempThetaWeak(:,j)' flip(tempThetaStrong(:,j)')],[tempBetaWeak(:,j)' flip(tempBetaStrong(:,j)')] );
     xlabel('Deflection Angle \theta');ylabel('Wave Angle \beta')
     title('\theta - \beta - \rmM \bfDiagram')
     % Mach Key
     text(45.4,70,'20','FontSize',7.5)
     text(41.3,67.3,'5','FontSize',7.5)
     text(34.2,65,'3','FontSize',7.5)
     text(23.2,63.5,'2','FontSize',7.5)
     text(12.2,65,'1.5','FontSize',7.5)
     text(5.3,69,'1.25','FontSize',7)
     text(.55,77.5,'1.05','FontSize',7)
     hold on
end


%% Problem 2

% Initialize AoA, half angles, Mach # 
alpha = 5;   % deg
epsilon1 = 10;  % deg
epsilon2 = 5;   % deg
Mach = 2;

% Call DiamondAirfoil.m
[c_l,c_dw] = DiamondAirfoil(Mach, alpha, epsilon1, epsilon2);

% Print lift and wave drag coefficients
fprintf('----- Problem 2 -----')
fprintf('\n')
fprintf('\n')
fprintf('Sectional Lift coefficient C_l: %0.4f \n',c_l)
fprintf('\n')
fprintf('Sectional Wave Drag coefficient C_dw: %0.4f \n',c_dw)
fprintf('\n')



%% Problem 3 

fprintf('----- Problem 3 -----')
fprintf('\n')
fprintf('\n')
fprintf('(Figure 2, Figure 3)')
fprintf('\n')
fprintf('\n')
fprintf('\n')
fprintf('{NOTE: For Figures 2 and 3, Shock Expansion Theory results are represented  \n') 
fprintf(' with a solid line while Linearized Theory results are represented with a dashed line.}')  
fprintf('\n')

% Initialize alpha vector
alphaVec = 0:9;
% Mach vector
MVec = [2 3 4 5];
% Cl and Cdw vectors
c_lVec = zeros(length(alphaVec),length(MVec));
c_dwVec = zeros(length(alphaVec),length(MVec));
c_lLinVec = zeros(length(alphaVec),length(MVec));
c_dwLinVec = zeros(length(alphaVec),length(MVec));

% color vector
colorVec = ['r' 'b' 'c' 'm'];

% Loop thru DiamondAirfoil.m for M = 2,3,4,5 and Alpha vector 
for k = 1:length(MVec)
    
    for m = 1:length(alphaVec)
    [c_lVec(m,k),c_dwVec(m,k)] = DiamondAirfoil(MVec(k), alphaVec(m), epsilon1, epsilon2);
    % From Linearized Theory
    c_lLinVec(m,k) = (4*deg2rad(alphaVec(m))) / sqrt(MVec(k)^2 - 1) ;
    c_dwLinVec(m,k) = (4*deg2rad(alphaVec(m))^2) / sqrt(MVec(k)^2 - 1) ;
    end
    
    figure(2);
    plot(alphaVec,c_lVec(:,k),colorVec(k),'LineWidth',2)
    hold on
    plot(alphaVec,c_lLinVec(:,k),colorVec(k),'LineStyle','--','LineWidth',2)
    ylim([0 0.5])
    xlabel('\alpha');ylabel('C_{l}')
    title('Coefficient of Lift vs. Angle of Attack')
    hold on
    
    figure(3);
    plot(alphaVec,c_dwVec(:,k),colorVec(k),'LineWidth',2)
    hold on 
    plot(alphaVec,c_dwLinVec(:,k),colorVec(k),'LineStyle','--','LineWidth',2)
    xlabel('\alpha');ylabel('C_{dw}')
    title('Coefficient of Wave Drag vs. Angle of Attack')
   
    
    
end

% Plot Cl vs Alpha
figure(2);
legend('M = 2','','M = 3','','M = 4','','M = 5','');
% Plot Cdw vs Alpha
figure(3);
legend('M = 2','','M = 3','','M = 4','','M = 5','');


% End


