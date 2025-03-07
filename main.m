% Author(s): Ashley Knight
% Assignment title: Project 2
% Purpose: Bottle rocket calculations and trajectory
% Creation date: 12/3/23
% Revisions:

close all; clear; clc;
%% Setup

const = getConst();
Thrust = 0;


tspan = [0 5];
xh = cosd(const.theta0);
zh = sind(const.theta0);
VAir_i = const.VB - const.VWater_i;
PTotal = const.Patm + const.PGage_i;
mAir_i = (PTotal * VAir_i) / (const.R * const.TAir_i);
Vx = const.v0 * cosd(const.theta0);
Vz = const.v0 * sind(const.theta0);



mRocket_i = const.mB + (const.rhoWater * const.VWater_i) + mAir_i;
initialConditions = [const.x0 const.z0 Vx Vz mAir_i VAir_i mRocket_i]';
initialConditions = double(initialConditions);

[t, X] = ode45(@(t, y) g_fun(t, y, const), tspan, initialConditions);
[output] = g_fun(t, X, const);

for i = 1:length(X(:,1))
    if X(i,2) < 0
       X = X(1:(i-1),:);
        t = t(1:i-1);
        break
    end
end


fprintf('The max height is %3.3d\n', double(max(X(:,2))))
fprintf('The max distance is %3.3d\n', double(max(X(:,1))))

thrustVec = createThrustVector(X,const);

figure(1)
plot(X(:,1), X(:,2), 'r')
xlabel('Distance [m]')
ylabel('Height [m]')
title('Launch Trajectory')


figure(2)
hold on
plot(t, thrustVec, 'r')
plot(verification.time, verification.thrust, 'b')
xlabel('Time [s]')
ylabel('Thrust [N]')
title('Thrust vs Time')
xlim([0 0.5])
hold off







function [thrustVec] = createThrustVector(X,const)
 
VAir_i = const.VB - const.VWater_i;
PTotal = const.Patm + const.PGage_i;
mAir_i = (PTotal * VAir_i) / (const.R * const.TAir_i);
  
  
   for i = 1:length(X(:,1))
  
       mAir = X(i,5);
       VAir = X(i,6);
  
       if VAir < const.VB
           % Water thrustVec Phase
           P_bottle = PTotal * (VAir_i / VAir) ^ const.gamma;
           thrustVec(i) = 2 * const.Cdis * const.aThroat * (P_bottle - const.Patm);
           
       else
           % All Water is gone
           PEnd = PTotal * (VAir_i / const.VB) ^ const.gamma;
           PAir = PEnd * (mAir / mAir_i) ^ const.gamma;
           % Test Case: Patm - PAir
           rhoAir = mAir / const.VB;
           TAir = PAir / (rhoAir * const.R);
           if PAir > const.Patm
               % Pressurized Air thrustVec(i)
               PCr = PAir * (2 / (const.gamma + 1)) ^ (const.gamma / (const.gamma - 1));
                   if PCr > const.Patm
                       % If p* > pa, the flow is choked (exit Mach number Me = 1)
                       Pe = PCr;
                       Te = (2 / (const.gamma + 1)) * TAir;
                       rho_e = Pe / (const.R * Te);
                       ve = sqrt(const.gamma * const.R * Te);
                   else
                       % If p* < pa, the flow is not choked
                       Pe = const.Patm;
                       Me = sqrt((2 / (const.gamma - 1)) * (((PAir / const.Patm) ^ ((const.gamma - 1) / const.gamma)) - 1));
                       Te = TAir / (1 + ((const.gamma - 1) / 2) * Me ^ 2);
                       rho_e = Pe / (const.R * Te);
                       ve = Me * sqrt(const.gamma * const.R * Te);
                   end
               mflow_Air = -1* (const.Cdis * rho_e * const.aThroat * ve);
               thrustVec(i) = (-1 * mflow_Air * ve) + (Pe - const.Patm) * const.aThroat;
               %phaseChange(2) = i;
              
           else
               % Balistic Phase -> Only gravity acting on rocket
               thrustVec(i) = 0;
           end
       end
   end
end

function [results] = g_fun(t,State,const)
   %Project2Equations defines the differential equations for the flight of
   % the bottle rocket
  Patm = 12.1 * 6894.75729;
   VAir_i = const.VB - const.VWater_i;
PTotal = const.Patm + const.PGage_i;
mAir_i = (PTotal * VAir_i) / (const.R * const.TAir_i);
   %% Read in initial values
   xPos = State(1);
   zPos = State(2);
   Vx = State(3);
  Vz = State(4);
   mAir = State(5);
   VAir = State(6);
   mRocket = State(7);
  
   vnorm = sqrt(Vx ^ 2 + Vz ^ 2);
   xh = cosd(const.theta0);
   zh = sind(const.theta0);
  
   
  
   Drag = (const.rhoAir / 2) * (vnorm ^ 2) * const.CD * const.aBottle;
   %% Thermo
  
   if VAir < const.VB
       % Water Thrust Phase
       P_bottle = PTotal * (VAir_i / VAir) ^ const.gamma;
       Thrust = 2 * const.Cdis * const.aThroat * (P_bottle - Patm);
       VDotAir = const.Cdis * const.aThroat * sqrt((2 / const.rhoWater) * (PTotal * ((VAir_i / VAir) ^ const.gamma) - Patm));
       mflow_r = - const.Cdis * const.aThroat * sqrt(2 * const.rhoWater * (P_bottle - Patm));
       mflow_Air = 0;
  
   else
       % All Water is gone
       PEnd = PTotal * (VAir_i / const.VB) ^ const.gamma;
       PAir = PEnd * (mAir / mAir_i) ^ const.gamma;
       % Test Case: Patm - PAir
       rhoAir = mAir / const.VB;
       TAir = PAir / (rhoAir * const.R);
       if PAir > Patm
           % Pressurized Air Thrust
           PCr = PAir * (2 / (const.gamma + 1)) ^ (const.gamma / (const.gamma - 1));
              
               if PCr > Patm
                   % If p* > pa, the flow is choked (exit Mach number Me = 1)
                   Pe = PCr;
                   Te = (2 / (const.gamma + 1)) * TAir;
                   rho_e = Pe / (const.R * Te);
                   ve = sqrt(const.gamma * const.R * Te);
               else
                   % If p* < pa, the flow is not choked
                   Pe = Patm;
                   Me = sqrt((2 / (const.gamma - 1)) * (((PAir / Patm) ^ ((const.gamma - 1) / const.gamma)) - 1));
                   Te = TAir / (1 + ((const.gamma - 1) / 2) * Me ^ 2);
                   rho_e = Pe / (const.R * Te);
                   ve = Me * sqrt(const.gamma * const.R * Te);
                  
               end
           mflow_Air = -1* (const.Cdis * rho_e * const.aThroat * ve);
           Thrust = (-1 * mflow_Air * ve) + (Pe - Patm) * const.aThroat;
           mflow_r = mflow_Air;
           VDotAir = 0; % -> No Air Flowing out
       else
           % Balistic Phase -> Only gravity acting on rocket
           Thrust = 0;
           VDotAir = 0;
           mflow_r = 0;
           mflow_Air = 0;
       end
   end
  
   
   xDot = Vx;
   zDot = Vz;
      
   if sqrt((xPos - const.x0) ^ 2 + (zPos - const.z0) ^ 2) > const.ls
       % If Rocket is still on stand, heading does not change
       xh = xDot / vnorm;
       zh = zDot / vnorm;
   end
  
   xVelDot = (Thrust - Drag) * (xh / mRocket);
   zVelDot = (Thrust - Drag) * (zh / mRocket) - const.g; 
  
   results = [xDot zDot xVelDot zVelDot mflow_Air VDotAir mflow_r]';
end
function const = getConst()
  
   const.g = 9.8;                                % Gravity [m/s^2]
   const.Cdis = 0.8;                               % discharge coefficient
   const.rhoAir = 0.961;                      % Ambient air density [kg/m^3]
   const.ls = 0.5;                               % length of test stand [m]
   const.tspan = [0 5];                      % time [s]
   VB = 0.002;                        % volume of empty bottle [m^3]
   Patm = 12.1;
   const.Patm = 12.1 * 6894.75729;               % atmospheric pressure [pa]
   const.gamma = 1.4;                            % ratio of specific heats for air
   const.rhoWater = 1000;                        % density of water [kg/m^3]
   dThroat = 0.021; 
                                                   % diameter of throat [m]
                                              
  dBottle = 0.105;                        % diameter of bottle [m]
   const.aThroat = pi * (dThroat / 2)^2;         % Area of Throat [m^2]
   const.aBottle = pi * (dBottle / 2)^2;         % Area of Bottle [m^2]
   const.dThroat =dThroat;  
   const.dBottle = dBottle;
   const.R = 287;                                % gas constant of air [J/kgK]
   const.mB = 0.15;                         % mass of empty 2-liter bottle with cone and fins [kg]
   const.TAir_i = 300;                      % initial temperature of air [K]
   const.x0 = 0.0;                         % initial horizontal distance [m]
   const.z0 = 0.25;                        % initial vertical height [m]
  
   const.v0 = 0.0;                         % initial velocity of rocket [m/s]
   const.theta0 = 42;                      % initial angle of rocket
    
   const.PGage_i = (52+12.1) * 6894.75729;        % initial gage pressure of air in bottle [pa]
   VWater_i = 0.001;                  % initial volume of water inside bottle [m^3]
   const.CD = 0.5;  % drag coefficient
   const.VAir_i = VB - VWater_i;
   const.VB = VB;
   const.VWater_i = VWater_i;
end
