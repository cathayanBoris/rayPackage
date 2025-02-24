function sig = sigmafit(kin,lin,Ei)
% function sig = sigmafit(kin,lin,Ei)
% calculate TRW frequency, sig, from wavevector components (k,l) 
%     for environmental parameters Ei
%
% rearrange sigma = NB*KXgh/(Kb*tanh(zeta)); where zeta includes sigma
% times K0/K0:
% sigma/(NB*KXgh/K0) = K0/(Kb* tanh(zeta)); 
% sig1= NB*KXgh/K0;         % define sig1 =freq for Kb=K0 & tanh(mu*h)=1
%
% sigma/sig1 = K0/(Kb* tanh(zeta));   % (both sides unitless)
%
% should work for any coordinate orientation 
% KXgh = hy*k - hx*l;       % 'K cross grad h'
% K0=sqrt(k^2 + l^2); 
% Kbsq= K0^2 + beta*k/sigma; %could test if Kbsq>=0, to exclude imaginary Kb
% Kb = sqrt(Kbsq);
% zeta= (NB/f)*Kb*h;
%
% define ratio  rat= sig1/sigma;  % rat= 0.01 :0.0001 :1.25
% rat=0.05: 0.000001: 1.25; % fine comb 'ratio', inverse for best res near 1
% find best fit lhs vs rhs,  for sigtest= sig1/rat;
%

h =  Ei(1);   hx=Ei(2);   hy=Ei(3);  
% hxx =Ei(4);  hxy=Ei(5);  hyy=Ei(6);  % not needed in sigmatest
beta=Ei(9);    f=Ei(8);    NB=Ei(7); 

if(length(kin) ~= length(lin))
    'input k and l vectors should have same length'
    return
end %if

sig=zeros(length(kin),1);  %preallocate space (column vector)

for nk=1:length(kin)
    k=kin(nk);
    l=lin(nk);   %loop one pair of (k,l) at a time


KXgh = hy*k - hx*l;          % 'K cross grad h', one (k,l) at a time
K0=sqrt(k*k + l*l);
sig1 = NB*KXgh/K0;            % define sig1 (1/s) for (k,l), when tanh=1

%%%%%%%%%%%%%%%%%%%%%%%
rtest = 0.01E-5; % choose rtest here; % 0.02E-5 runs tolerably quickly
   %1E-5 runs quickly 
   %0.005E-5 runs very slow (~10 minutes for 2x3d paths)
rat=0.1: rtest: 1.25; % 'ratio', use inverse for best resolution near 1
% rat=0.1: 0.000001: 1.25; % very fine comb 'ratio', inverse for best res near 1
% sigtest= sig1.*(1./rat);   % array depends inversely on rat
    % need to test many sigs between ~(0.8 to 10)*(NB*KXgh/K0), units (1/s)
Lr=length(rat); 
%%%%%%%%%%%%%%%%%%%%%%%

lhs=  zeros(Lr,1);  %initialize array sizes (column vectors)
rhs=  zeros(Lr,1);
% diff= zeros(Lr,1);   % calculated just once after the for jr loop

for jr=1:Lr 
%     sigr= sig1.* (1/rat(jr)); % units (figure1/s) one value of ratio at a time
    sigr= sig1/rat(jr);         % units (1/s) scalar, not a vector array
    Kbr = sqrt(K0*K0 + (beta*k)/sigr); % 'r' reminds sigr dependence
    zetar = (h*NB/f)*Kbr;               % 'r' reminds sigr dependence

    lhs(jr,1) = sigr/sig1;   % steps thru lhs=(sig1*1/rat)*sig1 = 1/rat
       % PROBLEM when KXgh->0, no restoring force, sig1->0, lhs->inf
          % and rhs= (K0/Kb)*(1/tanh(zeta), that depend on sigr
    rhs(jr,1) = K0/(Kbr* tanh(zetar)); 

end %for jr

% diff    = lhs-rhs; 
% absdiff = abs(diff);
% [~, index] = min(absdiff); 
[~, index] = min(abs(lhs-rhs));  
sig(nk,1) = lhs(index).*sig1;   %the best fit is at this index, lhs=1/rat
  
end %for nk loop

% % test plot rhs and lhs versus jr=1:Lr
% figure(9)
% jr=1:Lr;
% plot(jr,lhs,'k.', jr,rhs,'r.');
% title('lhs=1/rat (blk) & rhs=K0./(Kbr.* tanh(zetar)) (red)')
% xlabel('index for 1/rat')
% ylabel('lhs and rhs (dimensionless)')





