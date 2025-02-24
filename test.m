clear all
% load UGOS1sites.mat
% figure
% scatter(cpies.lon,cpies.lat)
%%
% a=find(cpies.lon < -82.25);
% b = length(a);
% lonInput(1:b) = cpies.lon(a);
% lonInput(b+1) = -88.006;
% latInput(1:b) = cpies.lat(a);
% latInput(b+1) = 26.2385;
% cpies.site(b+1) = "C02";

% Stones -90.7920, 26.4040
% I2 -89.9711, 27.2281
% C02 -88.006 26.2385
% Another I2 -90 27.2


% plate center
% -90.0083, 27.4917

% Instrument Locations
% lonList = [-89.2 -89.2 -89.165 -89.2 -89.35 -89.5 -89.65 -89.8 -89.98 -90.5];
% latList = [27.502 27.361 27.22 27.08 26.902 26.719 26.542 26.36 27.23 26.75];
% 
% 25pts
% lonList2 = [-90.9 -90.85 -90.8 -90.65 -90.6 -90.5 -90.45 -90.35 -90.2 -90.1 -89.95 ...
%     -89.8 -89.7 -89.55 -89.45 -89.4 -89.3 -89.2 -89.15 -89.1 -89.05 -89 -88.95 -88.9 -88.85];
% latList2 = [26.2 26.3 26.4 26.5 26.65 26.75 26.9 27.05 27.15 27.2 27.2 27.3 27.35 27.4 27.45 ...
%     27.55 27.5 27.45 27.35 27.3 27.2 27.15 27.05 26.95 26.9];

lonInput = [-89.15 -89.1 -89.05];
latInput = [26.75:-0.25:26.25];

needBathy = 0;
smoothingParameterInKM = 20;
initialPeriodInDays = [10:10:40];
initialWavelengthInKM = [90 125 200];
numberOfDaysToRun = 15;
timestepInHr = 1/8;
% forward 1 or backward -1
forwardBackward = 1;
%%
for ii = 1:length(lonInput)
    currLon = lonInput(ii);
    currLat = latInput(ii);
    for jj = 1:length(initialPeriodInDays)
        for kk = 1:length(initialWavelengthInKM)
            if needBathy == 0 && ii == 1 && jj == 1 && kk == 1
                needBathy = 1;
            else
                needBathy = 0;
            end
            % TRWrunvars(lonInput(ii),latInput(ii),initialPeriodInDays(jj),InitialWavelengthInKM(kk),numberOfDaysToRun,timestepInHr,forwardBackward,needBathy,smoothingParameterInKM);
            runTRWpath2ways(lonInput(ii),latInput(ii),initialPeriodInDays(jj),initialWavelengthInKM(kk),numberOfDaysToRun,timestepInHr,forwardBackward,needBathy,smoothingParameterInKM) 

        end
    end
end
%%
%        _                        
%        \`*-.                    
%         )  _`-.                 
%        .  : `. .                
%        : _   '  \               
%        ; *` _.   `*-._          
%        `-.-'          `-.       
%          ;       `       `.     
%          :.       .        \    
%          . \  .   :   .-'   .   
%          '  `+.;  ;  '      :   
%          :  '  |    ;       ;-. 
%          ; '   : :`-:     _.`* ;
%       .*' /  .*' ; .*`- +'  `*' 
%       `*-*   `*-*  `*-*'