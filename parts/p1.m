%% Fundamentals of GPS - Lab 3 - Part 1: Improved GPS Positioning

clear
clc
close all

%% Part A

% Import Trimble & Ionosphere Correction Data
load('RCVRT.mat')
load('RCVRT_data.mat')
load('IONO_CORR.mat')

% Antenna Truth 
ecefTrue = [423203.359 -5361678.541 3417280.681];
llaTrue = ecef2lla(ecefTrue);

figure
geoplot(llaTrue(:,1),llaTrue(:,2),'*')

% Instantiate Receiver
rcvr = gnssReceiver();

% Calculate GPS Position
RCVRTSteps = length(RCVRT);

% Spheroid Model
wgs84 = wgs84Ellipsoid('meter');

% Log Preallocation
ecefT = zeros(RCVRTSteps,3);
llaT = zeros(RCVRTSteps,3);
az = NaN(RCVRTSteps,10);
el = NaN(RCVRTSteps,10);

for i = 1:RCVRTSteps

    estT = rcvr.pv3D(RCVRT{i}.L1.psr, RCVRT{i}.L1.dopp, ...
        RCVRT{i}.L1.svPos, RCVRT{i}.L1.svVel, RCVRT{i}.L1.clkCorr, 1);

    % Log ECEF Position
    ecefT(i,:) = estT.pos;

        % Convert to LLA
    estT.lla = ecef2lla(estT.pos');

    % Log LLA Position    
    llaT(i,:) = estT.lla;

    % AZ & EL Calculations
    [RCVRT{i}.L1.az, RCVRT{i}.L1.el, ~] = ecef2aer(RCVRT{i}.L1.svPos(:,1), ...
        RCVRT{i}.L1.svPos(:,2), RCVRT{i}.L1.svPos(:,3), ...
        estT.lla(1), estT.lla(2), estT.lla(3), wgs84); %#ok<SAGROW> 

end

hold on
geoplot(llaT(:,1),llaT(:,2),'*')

%% Part B

% Window Size
M = 4200; % 7 min.

% IF Preallocation
IF = cell(RCVRTSteps,1);

% Initialization
IF{1}.S.psr = RCVRT{1}.L1.psr;
% 
% Log Preallocation
llaS = zeros(RCVRTSteps-1,3);

for i = 1:RCVRTSteps-1
    numSV = length(RCVRT{i+1}.L1.SVs);
    numSV_ = length(RCVRT{i}.L1.SVs);

    for j = 1:numSV

        psr = RCVRT{i+1}.L1.psr(j);
        carr = RCVRT{i+1}.L1.carr(j);
        
        if numSV > numSV_
            psr_ = RCVRT{i+1}.L1.psr(j);
            carr_ = RCVRT{i+1}.L1.carr(j);
        else
            psr_ = IF{i}.S.psr(j);
            carr_ = RCVRT{i}.L1.carr(j);
        end

        IF{i+1}.S.psr(j) = (1/M)*psr +  ((M-1)/M)*(psr_ + (carr-carr_));

    end

    estST = rcvr.pv3D(IF{i}.S.psr, RCVRT{i}.L1.dopp, ...
        RCVRT{i}.L1.svPos, RCVRT{i}.L1.svVel, ...
        RCVRT{i}.L1.clkCorr, 2);

    % Convert to LLA
    estST.lla = ecef2lla(estST.pos');

    % Log LLA Position    
    llaS(i,:) = estST.lla;

end

hold on
geoplot(llaS(:,1),llaS(:,2),'*')

%% Part C

alpha = [iono_corr.alpha_0; iono_corr.alpha_1; iono_corr.alpha_2; ...
    iono_corr.alpha_3];
beta = [iono_corr.beta_0; iono_corr.beta_1; iono_corr.beta_2; ...
    iono_corr.beta_3];

for i = 1:RCVRTSteps
    numSV = length(RCVRT{i}.L1.SVs);
    RCVRT{i}.L1.IEMD = zeros(numSV,1); %#ok<SAGROW>

    GPStime = RCVRT{i}.L1.gpsTime;
    Latitude = llaT(i,1);
    Longitude = llaT(i,2);
    

    for j = 1:numSV
        Azimut = RCVRT{i}.L1.az(j);
        Elevation = RCVRT{i}.L1.el(j);

        [~, RCVRT{i}.L1.IEMD(j), ~, ~] = ... 
        IonosphericDelay(GPStime,Latitude,Longitude,Azimut,Elevation,alpha,beta);
    end

end

% Log Preallocation
ecefIEM = zeros(RCVRTSteps,3);
llaIEM = zeros(RCVRTSteps,3);

for i = 1:RCVRTSteps

    psr = RCVRT{i}.L1.psr - RCVRT{i}.L1.IEMD;

    estIEMT = rcvr.pv3D(psr, RCVRT{i}.L1.dopp, ...
        RCVRT{i}.L1.svPos, RCVRT{i}.L1.svVel, RCVRT{i}.L1.clkCorr, 1);

    % Log ECEF Position
    ecefIEM(i,:) = estIEMT.pos;

        % Convert to LLA
    estIEMT.lla = ecef2lla(estIEMT.pos');

    % Log LLA Position    
    llaIEM(i,:) = estIEMT.lla;

end

hold on
geoplot(llaIEM(:,1),llaIEM(:,2),'*')

%% Part D
% NOTE: The IF pseudorange becomes incredibly noisy, therefore, the
% estimate is skewed.

% L-Band Frequencies
L1 = 1575.42e6;
L2 = 1227.60e6;
L3 = 1381.05e6;
L5 = 1176.45e6;

% Log Preallocation
llaDF = zeros(RCVRTSteps,3);

for i = 1:RCVRTSteps

    % Calculate Dual Frequency IF
    likeSV = ismember(RCVRT{i}.L1.SVs, RCVRT{i}.L2.SVs);
    likeSVidx = find(likeSV);

    IF{i}.DF.psr = ( ( RCVRT{i}.L1.psr(likeSVidx) * (L1^2/(L1^2 - L2^2)) ) ...
        - ( RCVRT{i}.L2.psr * (L2^2/(L1^2 - L2^2)) ) );

    estDFT = rcvr.pv3D(IF{i}.DF.psr, RCVRT{i}.L2.dopp, ...
        RCVRT{i}.L1.svPos(likeSVidx,:), RCVRT{i}.L1.svVel(likeSVidx,:), ...
        RCVRT{i}.L2.clkCorr, 2);

    % Convert to LLA
    estDFT.lla = ecef2lla(estDFT.pos');

    % Log LLA Position    
    llaDF(i,:) = estDFT.lla;

end

hold on
geoplot(llaDF(:,1),llaDF(:,2),'*')
legend('Truth','Original','Carrier Smoothed','IEM IF','Dual Frequency IF','Location','northwest')
title('IF & Smoothed Pseudorange Solution Comparison')

geobasemap satellite

% Rearrange Plot Order
plot = get(gca,'Children');
order = [3 2 4 1 5];
plot = plot(order);
set(gca,'Children',plot)
