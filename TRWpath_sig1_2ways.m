function PVAR =TRWpath_sig1_2ways(initialPeriodInS ,initialWavelengthInM ,xq,yq, Ei, dt_hr,nsteps,...
    Bgxy, Bgrid, smoothPara)

global meandx meandy originLon originLat omega

% function PVAR =TRWpath(T0,LAM0, xq,yq, Ei, dt_hr,nsteps,Bath)
%
% calculates TRW path-including beta;
% NB depends on bottom depth h,
% added NB(h) dependence, with NB=const above h=bottom depth
%    easy, done within each step of Runge-Kutta loop
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ei is an INPUT to TRWpath,
% AND Ei also gets updated 2X, at mid-point and end-point for each R-K step

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% will output calculated PVAR all along path, where
% PVAR=[xp,yp,kp,lp,K0p,LAMkmp,sigp,Tdayp,cgxp,cgyp,hp,hxp,hyp];
%

% initial inputs needed, but not used on return trip
%         T0day =  TRW period (days)
sig0=2*pi/(initialPeriodInS );   % units (1/s), exact from T0 (sec)

%        LAM0 = TRW wavelength (m)
K0=2*pi/(initialWavelengthInM );  % units (1/m), from LAM0 (m)


dt_s= dt_hr*3600;    % units =(sec)  % dt_hr is now an input var
% nsteps ; % # steps to track, one-way = input variable calculated outside this function

X0= xq;  Y0= yq;  %starting position (m) was specified in "running_path"

% initial Ei values are input from running_path.m
h = Ei(1); hx = Ei(2); hy = Ei(3);
hxx = Ei(4); hxy = Ei(5); hyy = Ei(6);
NB = Ei(7); f = Ei(8); beta = Ei(9);

iUPDN=-1;  % normal case (-1) is DNSLOPE phase propagation and cg upslope

%%%%%% will need whole gridded Bgrid before entering R-K loop %%%%%%%%%%
%%%%%%%%%%%%%%% whole gridded fields of bathymetry & gradients %%%%%%%%%
% ifplots=0;
% [Bgxy, Bgrid] =bathy_sim(Bath,ifplots); %instead, these are input variables
% to avoid recalculating fine-scale bathy
% Bgxy are the x and y grid axes, stored in 1 long row,
%       then y gets transposed into a column vector
% Bgrid has all six gridded fields

% offload xyh combined grids
[~,ly]=size(Bgxy); % = 2lyg
xg=Bgxy(:,1:ly/2);
yg=Bgxy(:,ly/2+1:ly)';

[~,lx]=size(Bgrid); % 6lxg
lx = lx/6;

hg =Bgrid(:,1:lx)';
hxg=Bgrid(:,lx+1:2*lx)';
hyg=Bgrid(:,2*lx+1:3*lx)';
hxxg=Bgrid(:,3*lx+1:4*lx)';
hxyg=Bgrid(:,4*lx+1:5*lx)';
hyyg=Bgrid(:,5*lx+1:6*lx)';

%%%%%%%%%%%%%%%%%%% calc starting wave vector, (k,l)  %%%%%%%%%%%%%%%%%
% use 10X finer-comb fit for 1st point on path,
% append "p" in klcirclefitp() for this slow one-time fit

[k0,l0,sigFitted]=klcirclefitp(sig0,K0,Ei,iUPDN); %function name has appended "p" for path

% redefine sig0 so it fits (k,l) even for long waves, small K0
% DON'T DO THE FOLLOWING ANY MORE
%  sig2=sigmafit(k0,l0,Ei);    % sig2 fit is (1/5M=0.2ppm) in sigmafit
%  % already output at this point, sig1 fit is finer (1/128M) in klcirclefitp
% keyboard    %testing point
%[k0,l0, sig1]=klcirclefitp(sig0,K0,Ei,iUPDN); sigtest=sigmafit(k0, l0, Ei);
%  (sig1-sig0)/sig0, (sigtest-sig0)/sig0
% gives -8.1918e-08,   -1.1830e-07  ... sig1 is 35% closer to sig0
% so I changed sigmafit output sig1 to use rhs instead of lhs that had
% been identically set =sig0

% DRW verified that sigtest=sigmafit(k0,l0,Ei) agrees with sig1 within 0.1ppm
% simple sanity check on sig1
dsfrac=(sigFitted-sig0)/sig0;
if(abs(dsfrac)>1E-5)  % passes 1E-5, try 1E-6
    'TRWpath initialization mismatch sig,k,l, SKIP'
    PVAR = 999;
    return
    % keyboard  % I suspect there would be an error if you get here
else
    %%%%  NOTE change to EXACT k,l,sig1, Ei fit from rhs of klcirclefitp %%%
    sig0=sigFitted; % for normal initialization, make this tiny <1ppm adjustment
    T0day=2*pi/sigFitted/86400;
    % testing forward+reverse path this initialization choice makes no difference
end %if

% k,l,sig0 should now exactly fit the dispersion curve because changing
% to sig0=sig1 output by klcirclefitp rhs is exact fit output for k,l,Ei

%%%%%%%%%% know starting wave vector and freq, (k,l,sig0=sig1)  %%%%%%%%%%
% [X0,Y0] % temporary test-value checks
% [k0, l0]Bath
% Kbsqr =(k0*k0 + l0*l0 + k0*beta/sig0); % this 1st value is printed just to test
%return

% npoints= 2*nsteps-1;    % initialize space also for reverse path
npoints = nsteps; % now it forbids the fw-bw tracing, only one directional 
% initialize variables - space for vars that we may want to plot along path
%   VECTOR LIST                     % SCALAR 1st-values LIST
xPath = zeros(npoints,1);               xPath(1,1)=X0;   % append 'p' for values along path
yPath = zeros(npoints,1);               yPath(1,1)=Y0;
kPath = zeros(npoints,1);               kPath(1,1)=k0;
lPath = zeros(npoints,1);               lPath(1,1)=l0;
K0Path= zeros(npoints,1);               K0Path(1,1)= K0;   %1st value might differ ??
periodPathinDays = zeros(npoints,1);    periodPathinDays(1,1)=T0day;
cgxPath = zeros(npoints,1);      % get 1st cgx, cgy vals in 1st pass of for-loop
cgyPath = zeros(npoints,1);
EiPath = zeros(npoints,9);              EiPath(1,:) = Ei(:,:);
% hPath = zeros(npoints,1);      hPath(1,:) = Ei(:,1);
% hxPath = zeros(npoints,1);      hxPath(1,:) = Ei(:,2);
% hyPath = zeros(npoints,1);      hyPath(1,:) = Ei(:,3);
% NBPath = zeros(npoints,1);      NBPath(1,:) = Ei(:,7);
% KbPath = zeros(npoints,1);     KbPath(1,1) = sqrt(kPath(1,1).^2 + lPath(1,1).^2 + kPath(1,1)*beta/sig0);


% wait until end to define PVAR:
% PVAR= [xp,yp,kp,lp,K0p,LAMkmpp,sigp,Tdayp,cgxp,cgyp,hp,hxp,hyp]; %wait until end

% % % the following subset are used only for initial approx sig sanity-check
% % Kb =sqrt(k0*k0 + l0*l0 + k0*beta/sig0);  % K with beta - uses input sig0
% % %kbs=k+ beta/(2*sig0) ;                % k with beta term - uses input sig0
% % KXgh = (hy*k0 - hx*l0);                 % "K cross grad h"
% % zeta= N*Kb*h/f;                   % argument in tanh and sinh
% % Tz = tanh(zeta);                  % the tanh term
% % sigapprox=N*KXgh/(Kb*Tz);
% % dsigrel= (sigapprox-sig0)/sig0;
% %

figure(100)
text(0, -1e5, ['dt is ', num2str(dt_hr,3)] )
text(0, -1.5e5, ['T0day=', num2str(T0day,3), '  λ0km=', num2str(initialWavelengthInM/1000,5)])
hold on


%%%%%%       %Runge-Kutta loop       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig=sig0;  %startoff value for R-K loop; after klcirclefit, set sig0=sig1
% keyboard   %uncomment for debugging

nreverse=nsteps+1; %this point is the halfway point
for ns=2:npoints    % note loop starts calc for 2nd point from 1st point
    disp([num2str((ns-1)),'/',num2str(npoints-1),' Complete']);
    if (ns==nreverse) %at halfway, reverse the path
        dt_s=-dt_s;
    end %if

    % x=xp(ns-1,1);    y=yp(ns-1,1);   % append 'p' for values along path
    k=kPath(ns-1,1);    l=lPath(ns-1,1);   % append 'p' for values along path

    %ready for 2nd order Runge-Kutta time-stepping method %Press,etal1987 p.550
    %   (for name of along-path step distance, just substitute s for Press's h
    %   Press's xn represent location variables (x,y)
    %   Press's yn represent function variables f
    % k1=sf'(xn,yn); %1st calc full step s=dt times the slope f' from init point
    % k2=sf'(xn+s/2, yn+(k1)/2);  %s=dt times the slope f' re-calc'd at the mid-point
    % y_{n+1} = (yn + k2) + O(s^3)
    %
    % the tentative mid-point vals, are NOT along final path
    % the k1 step uses postion and wavenumber evaluated at ns-1
    % the k2 step uses values halfway to next point at ns

    % the following four vars are equiv to Press's above set of k1*(1/2) values
    % Ei already is for step ns-1 location (x,y)
    % and (k,l,sig) are for step ns-1
    [cgx1,cgy1,dkdt1,dldt1] =TRWparms_point(k,l,sig,EiPath(ns-1,:)); %(k,l,sig,Ei) for step ns-1
    if(ns==2)  %1st step, save cgx and cgy
        cgxPath(1,1) =cgx1;
        cgyPath(1,1) =cgy1;
    end
    % x,y,k,l,sig get re-calculated at each time-step in this ns for-loop
    % vals at each R-K k1 point
    xk1Midpoint = xPath(ns-1,1)+ cgx1* dt_s/2;    % appended 'm' for 'mid-point'
    yk1Midpoint = yPath(ns-1,1)+ cgy1* dt_s/2;
    kk1Midpoint = kPath(ns-1,1)+ dkdt1*dt_s/2;
    lk1Midpoint = lPath(ns-1,1)+ dldt1*dt_s/2;
    %keyboard

    %%%%
    % also need all the gridded bathy parameters inside subroutine TRWpath
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % get environment parms at mid-point Ei= [h hx hy hxx hxy hyy NB f beta];
    %    NB at mid-point gets calculated from NBh(h);

    xg = xg(:,1);
    yg = yg(:,1);
    % examine the environmental parameters at the R-G midpoint

    latMidpoint = yk1Midpoint/60/meandy+originLat;
    fMidpoint =  2*omega*sind(latMidpoint);
    betaMidpoint = 2*omega*(sind(latMidpoint+0.01) - sind(latMidpoint-0.01))/(0.02*pi*6371000/180);
    EiMidpoint = Env_parms(xk1Midpoint,yk1Midpoint,xg,yg,hg',hxg',hyg',hxxg',hxyg',hyyg',fMidpoint,betaMidpoint); %from grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % NO LONGER do a new sigmafit at halfway point=midpoint, preparing for k2 step
    % sigm =sigmafit(kk1m,lk1m,Eim); %previously gave time-reversible paths better than sig0
    sigMidpoint= sig0;

    % the next four are equiv to Press's set of k2 values (xn+ s/2, yn+ k1/2)
    % use inputs (kk1m,lk1m,sigm,Ei)   all updated in k1 step
    [cgx2,cgy2,dkdt2,dldt2] =TRWparms_point(kk1Midpoint,lk1Midpoint,sigMidpoint,EiMidpoint);

    % the new set of x,y,k,l values along the path, using mid-point gradients
    xPath(ns,1)= xPath(ns-1,1) +  cgx2*dt_s;
    yPath(ns,1)= yPath(ns-1,1) +  cgy2*dt_s;
    kPath(ns,1)= kPath(ns-1,1) + dkdt2*dt_s;
    lPath(ns,1)= lPath(ns-1,1) + dldt2*dt_s;

    %preparing for next iteration step
    x=xPath(ns,1);
    y=yPath(ns,1);
    %%%% Ei calculated at next step (ns) of R-K integration %%%%%%%%%%%%%
    latFinal = y/60/meandy+originLat;
    fFinal =  2*omega*sind(latFinal);
    betaFinal = 2*omega*(sind(latFinal+0.01) - sind(latFinal-0.01))/(0.02*pi*6371000/180);
    EiPath(ns,:) = Env_parms(x,y,xg,yg,hg',hxg',hyg',hxxg',hxyg',hyyg',fFinal,betaFinal); %from grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k=kPath(ns,1);
    l=lPath(ns,1);
    % sigp(ns,1) =sigmafit(k,l,Ei);  % here we refit sig at next R-K point
    % sig=sigp(ns,1);                % to feed into next R-K step after ns;
    %  this allowed sig to drift, but it previously worked consistently in reversing dt paths
    sig=sig0;

    figure(100)
    if ns<nreverse
        plot(xPath(ns,1),yPath(ns,1),'b+')
    elseif ns>=nreverse
        plot(xPath(ns,1),yPath(ns,1),'m+')
    end

    % % title('TRW path')
    % text(-90E3, -85E3, ['dt is ', num2str(dt_hr,3)] )
    % text(-90E3, -90E3, ['T0day=', num2str(T0day,3), '  LAM0km=', num2str(LAM0/1000,5)])
    % hold on

    if ns == npoints
        plot(xPath(1,1),yPath(1,1),'g*')
    end
    % move the text to here because this 1st point is plotted just once

    % do NOT need to recalc Ei here, already calculated at ns point in R-K step
    % must finalize cgxp, cgyp at this ns time-step
    % now at new (xp,yp) loc, calc cgxp(ns,1),cgyp(ns,1), inputting new sig to TRWparms
    %
    [cgxPath(ns,1),cgyPath(ns,1),~,~]=TRWparms_point(k,l,sig,EiPath(ns,:)); %now must save cgxp,cgyp
    %  identical to  = TRWparms_point(kp(ns,1),lp(ns,1),sigp(ns,1),Ei);
    
    %%%%%%%%%%%% abort mission if error %%%%%%%%%%%
    
    if (cgxPath(ns,1) == 999)
        nValid = 1:ns-1;
        break
    end

    % sig cannot drift from the constant sig0=sig1 value from klcirclefitp
    % % ADVANTAGE: this fitted sig stays exactly on the dispersion surface

    % not saved: zeta, tanh

    K0Path(ns,1)=  sqrt(kPath(ns,1).*kPath(ns,1) + lPath(ns,1).*lPath(ns,1));

    Kbsqr =(kPath(ns,1).^2 + lPath(ns,1).^2 + kPath(ns,1)*beta/sig);    % K with beta     - uses input sig
    gH = sqrt(EiPath(ns,2).^2+EiPath(ns,3).^2);
    if Kbsqr <= 0
        'in TRWpath_sig1_2ways Kb^2 is < 0'
        nValid = 1:ns-1;
        break % abort simulation
    end

    if dt_s > 0 % forward tracing only
        if x > 6E4 && EiPath(ns,2) > -0.0031 && gH > 0.01294
            'in TRWpath_sig1_2ways arrived at MS Fan'
            nValid = 1:ns-1;
            break % abort simulation
        end
        if(EiPath(:,1) <= 1500)
            'in TRWpath_sig1_2ways -1500m limit has been reached'
            nValid = 1:ns-1;
            break % abort simulation
        end
        if x < 2E4
            if gH > 0.036
                'in TRWpath_sig1_2ways arrived at Sigsbee Escarpment'
                nValid = 1:ns-1;
                break % abort simulation
            end
        end
    end

    if dt_s < 0 % backward tracing only
        if x > 2E5
            'in TRWpath_sig1_2ways arrived at eastern edge of the domain'
            nValid = 1:ns-1;
            break % abort simulation
        end
        if y < -1.05E5
            'in TRWpath_sig1_2ways arrived at southern edge of the domain'
            nValid = 1:ns-1;
            break % abort simulation
        end
    end

    hVector = EiPath(ns,2) + EiPath(ns,3)*1i;
    KVector = k + l*1i;
    dot = real(hVector)*real(KVector)+imag(hVector)*imag(KVector);
    angleInRadian = acos(dot/(abs(hVector)*abs(KVector)));

    if(angleInRadian > pi/2) % not robust, can be improved
        'in TRWpath_sig1_2ways wavenumber vector entered wrong quadrant'
        nValid = 1:ns-1;
        break % abort simulation
    end

    % sanity check
    % [k0,l0, sigFitted]=klcirclefitp(sig0,K0,Ei,iUPDN);

    % if dsfrac > 1E-5
    %     2*pi/sig0/86400
    %     2*pi/frequencyPath/86400
    %     'SigmaPath mismatch sig,k,l, TERMINATED'
    %     nValid = 1:ns-1;
    %     break % abort simulation
    % end

    % NBPath(ns,1) = EiPath(ns,7);
    % hxPath(ns,1) = EiPath(ns,2);
    % hyPath(ns,1) = EiPath(ns,3);
    % hPath(ns,1) = EiPath(ns,1);
    %
    % Kbsqr =(K0Path(ns,1).^2 + kPath(ns,1)*beta/sig);    % K with beta     - uses input sig
    % KbPath(ns,1) = sqrt(Kbsqr);
    %
    % zeta(ns)=hPath(ns,1)*KbPath(ns,1)*NBPath(ns,1)/f;
    % KXgh(ns)  =  (hyPath(ns,1)*kPath(ns,1) - hxPath(ns,1)*lPath(ns,1));
    % sigCheck(ns) = NBPath(ns,1).*KXgh(ns)./KbPath(ns,1)./tanh(zeta(ns));
    % 2*pi/sigCheck(ns)/86400

    % seems OK?
    nValid = 1:ns;

end % for ns

%%%%%%       %END of Runge-Kutta loop       %%%%%%%%%%%%%%%%%%%%%

% figure(999)
% subplot(5,1,1)
% plot(NBPath)
% ylabel('N 1/s')
% subplot(5,1,2)
% plot(KXgh)
% ylabel('KX∇h')
% subplot(5,1,3)
% plot(KbPath)
% ylabel('Kb')
% subplot(5,1,4)
% plot(tanh(zeta))
% ylabel('tanh(ζ=μh)')
% subplot(5,1,5)
% plot(2*pi./sigCheck/86400)
% yline(25)
% ylim([-5 5]+initialPeriodInS/86400)
% ylabel('T day')


startT = initialPeriodInS/86400;
startLAM = initialWavelengthInM/1000;
startLon = xq/60/meandx+originLon;
startLat = yq/60/meandy+originLat;
pathLon = xPath(nValid,1)/60/meandx+originLon;
pathLat = yPath(nValid,1)/60/meandy+originLat;
pathK0 = K0Path(nValid,:);
pathK = kPath(nValid,1);
pathL = lPath(nValid,1);
% pathCgx = cgxPath(nValid,1);
% pathCgy = cgyPath(nValid,1);
% pathCg = sqrt(cgxPath(nValid,1).*cgxPath(nValid,1)+cgyPath(nValid,1).*cgyPath(nValid,1));
pathEi = EiPath(nValid,:);
smoothing = smoothPara/1000;
totalSteps = nValid(end);
days = abs((totalSteps-1)*dt_hr/24);
timeHr = 0:abs(dt_hr):(totalSteps-1)*abs(dt_hr);

% upon exiting R-K loop, now have vals for the full path
PVAR = 1;
% PVAR= [xPath,yPath,kPath,lPath,K0Path,wavelengthPathinKM,frequencyPath,periodPathinDays,cgxPath,cgyPath,hPath,hxPath,hyPath];

% distp=sqrt((xp-X0).*(xp-X0)+(yp-Y0).*(yp-Y0));  %distance from start pos

% Dmax=max(distp);  % distance (m), plotting limits that the TRWs travel
% hmax=4000;  hmin=500; %bathy range
% write warnings if...
if ((max(periodPathinDays)>64)|| (max(EiPath(:,1))>=3950) || (min(EiPath(:,1))<=600))
    'Tday >64 or h reached near bathy slope limits'
end %if

varName = append(num2str(initialPeriodInS/86400),'D',num2str(initialWavelengthInM/1000),'KM',num2str(startLon),'E', ...
    num2str(startLat),'N',num2str(dt_hr),'HdT',num2str(days),'day', ...
    num2str(smoothing),'KMsmooth');

%forward/backward save .mat
if dt_hr>0
    filePath = 'C:\Users\Yan_J\OneDrive\Documents\MATLAB\GoM\TRW\rayTracing\forward\';
elseif dt_hr<0
    filePath = 'C:\Users\Yan_J\OneDrive\Documents\MATLAB\GoM\TRW\rayTracing\backward\';
end

save([filePath varName '.mat'],"startT","startLAM","startLon","startLat","pathLon", ...
    "pathLat","dt_hr","pathK0","pathK","pathL", ...
    "pathEi","smoothing","totalSteps","days","timeHr","-mat");


%%%%%% saved variables list %%%%%%
% startT - initial period in days
% startLAM - initial wavelength in km
% startLon, startLat - initial location of tracing in deg
% pathLon, pathLat - location at each step in tracing in deg
% dt_hr, length of timestep in hr
% pathK0,pathK, pathL - wavenumbers at each step in tracing in 1/m
% pathEi - environmental parameters, h,hx,hy,hxx,hxy,hyy,N,f,beta
% smoothing - bathymetry smoothing parameter in km
% totalSteps - number of valid entries
% timeHr -time record of ray tracing, a virtual clock in hr



% figure(4); hold on;
% subplot(3,1,1); plot(distp,Tdayp); hold on;
% xlabel('distance (m) from start position'); ylabel('Tday (days)')
% axis([0 Dmax 0 64 ]);
% subplot(3,1,2); plot(distp,LAMkmp); hold on;
% xlabel('distance (m) from start position'); ylabel('LAMkm (km)')
% axis([0 Dmax 0 400]);
% subplot(3,1,3); plot(distp,hp); hold on;
% xlabel('distance (m) from start position'); ylabel('h depth (m)')
% axis([0 Dmax hmin hmax]);
%
% figure(5); hold on; % want to plot K0p, kp,lp
% subplot(3,1,1); plot(distp, K0p); hold on;
% xlabel('distance from start position'); ylabel('K0p (1/m)')
% %         Kmin=min(K0p); Kmax=max(K0p);
% Kmin=-3E-4; Kmax=3E-4;  % for wavelengths as short as ~21km.
% axis([0 Dmax Kmin Kmax]);
% subplot(3,1,2); plot(distp, kp); hold on
% xlabel('distance from start position'); ylabel('kp (1/m)')
% %         Kmin=min(kp); Kmax=max(kp);
% axis([0 Dmax Kmin Kmax]);
% subplot(3,1,3); plot(distp, lp); hold on
% xlabel('distance from start position'); ylabel('lp (1/m)')
% %         Kmin=min(lp); Kmax=max(lp);
% axis([0 Dmax Kmin Kmax]);
%
% figure(6); hold on;  % want to plot cgxp, cgyp, |cgp|
% subplot(3,1,1); plot(distp, cgxp); hold on
% xlabel('distance from start position'); ylabel('cgxp (m/s)')
% %         cmin=min(cgxp); cmax=max(cgxp);
% cmin=-3.0; cmax=+3.0;
% axis([0 Dmax cmin cmax]);
% subplot(3,1,2); plot(distp, cgyp); hold on
% xlabel('distance from start position'); ylabel('cgyp (m/s)')
% %         cmin=min(cgyp); cmax=max(cgyp);
% axis([0 Dmax cmin cmax]);  %use same limits as for cgx
% subplot(3,1,3); plot(distp, sqrt(cgxp.*cgxp+cgyp.*cgyp)); hold on
% xlabel('distance from start position'); ylabel('|cg| (m/s)')
% axis([0 Dmax 0 cmax]);
%
% figure(7); hold on;
% plot vert decay scale 1/mu;  mu~ (N/f)*K; 1/mu=f*(LAMkmp*1e3)/(2*pi*N)
% A=f/(2*pi);
% muinvp=(A*LAMkmp)./NBh(hp));

% zetabp=hp.*mubp;
% subplot(3,1,1); plot(distp,1e-3./mubp); hold on;   % units km
% xlabel('distance from start position'); ylabel('trapping scale 1/mub (km)')
% axis([0 Dmax 0 6 ]);
% % plot ratio Kb/K0 =sqrt(K0p.*K0p+kp.beta./sigp);
% subplot(3,1,2); plot(distp, Kbp./K0p ); hold on;
% xlabel('distance from start position'); ylabel('ratio Kb/K0')
% axis([0 Dmax 0 1.2]);
% subplot(3,1,3); plot(distp,1./tanh(zetabp)); hold on;
% xlabel('distance from start position'); ylabel('1/tanh(zeta)')
% axis([0 Dmax 0 20]);
%
% figure(8)
% lambda, Cg, N, alpha, 1/mub

% xLimit = [0 abs((nsteps-1)*dt_hr)];
% subplot(3,2,1)
% if dt_hr>0
%     plot(time,pathK0); hold on;
%     xlim(xLimit);
% elseif  dt_hr<0
%     plot(-time,pathK0); hold on;
%     xlim(fliplr(-xLimit));
%     set(gca,'XDir','reverse');
% end
% xlabel('time (h) from start position'); ylabel('LAMkm (km)')
% grid on
%
% subplot(3,2,2)
% if dt_hr>0
%     plot(time, pathCg); hold on
%     xlim(xLimit);
% elseif  dt_hr<0
%     plot(-time, pathCg); hold on
%     xlim(fliplr(-xLimit));
%     set(gca,'XDir','reverse');
% end
% xlabel('time (h) from start position'); ylabel('|Cg| (m/s)')
% grid on
%
% subplot(3,2,3)
% if dt_hr>0
%     plot(time, pathAlpha); hold on
%     xlim(xLimit);
% elseif  dt_hr<0
%     plot(-time, pathAlpha); hold on
%     xlim(fliplr(-xLimit));
%     set(gca,'XDir','reverse');
% end
% xlabel('time (h) from start position'); ylabel('Alpha (m/m)')
% grid on
%
% subplot(3,2,4)
% if dt_hr>0
%     plot(time, trappingScale); hold on
%     xlim(xLimit);
% elseif  dt_hr<0
%     plot(-time, trappingScale); hold on
%     xlim(fliplr(-xLimit));
%     set(gca,'XDir','reverse');
% end
% xlabel('time (h) from start position'); ylabel('trapping scale 1/mub (km)')
% grid on
%
% subplot(3,2,5)
% if dt_hr>0
%     plot(time,pathN); hold on; ylim([0 1.2*max(pathN)]);
%     xlim(xLimit);
% elseif  dt_hr<0
%     plot(-time,pathN); hold on; ylim([0 1.2*max(pathN)]);
%     xlim(fliplr(-xLimit));
%     set(gca,'XDir','reverse');
% end
% xlabel('time (h) from start position'); ylabel('N')
% grid on

