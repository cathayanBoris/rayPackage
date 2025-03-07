clear all

if exist('needBathy','var') && needBathy == 1
else
    needBathy = 0;
end

lonInput = [];
latInput = [];

smoothingParameterInKM = 20;
initialPeriodInDays = [];
initialWavelengthInKM = [];
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