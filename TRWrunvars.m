function [lonInput,latInput,initialPeriodInDays,InitialWavelengthInKM,numberOfDaysToRun,timestepInHr,forwardBackward,needBathy,smoothingParameterInKM] ...
    =TRWrunvars(lonInput,latInput,initialPeriodInDays,InitialWavelengthInKM,numberOfDaysToRun,timestepInHr,forwardBackward,needBathy,smoothingParameterInKM)
% function TRWSTARTUP=TRWrunvars
% these are the input variables needed for runTRWpath2way
%   function runTRWpath2way(X0km,Y0km,T0day,LAM0km,ndaysrun,dt_hr,sgn_dt,...
%       ifcalcbathy,Bath)

% X0km=-64.5*sqrt(2);  Y0km=0;  % positi on near h~1100m at western edge 

% Hamilton -89.4, 26.2
% Stones -90.7920, 26.4040
% L7 -86.97, 25.52
% PIES -88.6177, 25.6971
% C1 -88.6072, 26.2389
% D2 -88.0151, 25.6969
% D3 -87.4157, 25.6985
% A2 -87.987, 27.318
% -87.987, 27.2583
% -87.897, 27.2083
% B1 -88.6092, 26.7789
% B2 -87.9956 26.7771
% I2 -89.9711, 27.2281 
% C02 -88.006 26.2385
% W-Stones -90.883 26.404
% position near h~2500m at western edge

% % Bath=BATH; %global variable BATH not working, so do the following
% Bath= [34E3,3000,3800,60E3,2E3,0.50E3];    %units (m) % 0.5 km bathy grid,
%             % all values must be integer multiples of last number=dx

% TRWSTARTUP= [X0km,Y0km,T0day,LAM0km,ndaysrun,dt_hr,sgn_dt,ifcalcbathy];

runTRWpath2ways(lonInput,latInput,initialPeriodInDays,InitialWavelengthInKM,numberOfDaysToRun,timestepInHr,forwardBackward,needBathy,smoothingParameterInKM) 

%  _   _            _ 
% | | | | ___ _   _| |
% | |_| |/ _ \ | | | |
% |  _  |  __/ |_| |_|
% |_| |_|\___|\__, (_)
%             |___/   


