%% Fundamentals of GPS - Lab 3 - Part 1: Improved GPS Positioning

clear
clc

%% Part A

% Import Trimble Data
load('RCVRT.mat')

% Instantiate Receiver
rcvr = gnssReceiver();

% Calculate GPS Position
RCVRTSteps = length(RCVRT);

% Log Preallocation
llaT = zeros(RCVRTSteps,3);

for i = 1:RCVRTSteps

    estT = rcvr.pv3D(RCVRT{i}.L1.psr, RCVRT{i}.L1.dopp, ...
        RCVRT{i}.L1.svPos, RCVRT{i}.L1.svVel, RCVRT{i}.L1.clkCorr, 1);

    % Convert to LLA
    estT.lla = ecef2lla(estT.pos');

    % Log LLA Position    
    llaT(i,:) = estT.lla;

end

figure
geoplot(llaT(:,1),llaT(:,2),'*')

%% Part B

% Window Size
% win_sz = 4200; % 7 min.
% 
% IF Preallocation
% IF = cell(RCVRTSteps,1);
% 
% Initialization
% IF{1}.S.psr = RCVRT{1}.L1.psr;
% 
% Log Preallocation
% llaS = zeros(RCVRTSteps,3);
% 
% for i = 1:RCVRTSteps-1

%     lenk = length(RCVRT{i}.L1.SVs);
%     lenk1 = length(RCVRT{i+1}.L1.SVs);
% 
%     if lenk1 > lenk
%         likeSV = ismember(RCVRT{i+1}.L1.SVs, RCVRT{i}.L1.SVs);
%         likeSVidx = find(likeSV);
%     else
%         likeSV = ismember(RCVRT{i}.L1.SVs, RCVRT{i+1}.L1.SVs);
%         likeSVidx = find(likeSV);
%     end
% 
%     IF{i+1}.S.psr = ( 1/win_sz ) * RCVRT{i+1}.L1.psr(likeSVidx) + ( (win_sz-1)/win_sz ) ...
%         * ( IF{i}.S.psr + RCVRT{i+1}.L1.carr(likeSVidx) - RCVRT{i}.L1.carr );

%     estST = rcvr.pv3D(IF{i}.DF.psr, RCVRT{i}.L2C.dopp, ...
%         RCVRT{i}.L1.svPos(likeSVidx,:), RCVRT{i}.L1.svVel(likeSVidx,:), ...
%         RCVRT{i}.L2C.clkCorr, 2);
% 
%     % Convert to LLA
%     estST.lla = ecef2lla(estST.pos');
% 
%     % Log LLA Position    
%     llaS(i,:) = estST.lla;

% end


%% Part C

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
    likeSV = ismember(RCVRT{i}.L1.SVs, RCVRT{i}.L2C.SVs);
    likeSVidx = find(likeSV);

    IF{i}.DF.psr = ( ( RCVRT{i}.L1.psr(likeSVidx) .* (L1^2/(L1^2 - L2^2)) ) ...
        - ( RCVRT{i}.L2C.psr .* (L2^2/(L1^2 - L2^2)) ) );

    estDFT = rcvr.pv3D(IF{i}.DF.psr, RCVRT{i}.L2C.dopp, ...
        RCVRT{i}.L1.svPos(likeSVidx,:), RCVRT{i}.L1.svVel(likeSVidx,:), ...
        RCVRT{i}.L2C.clkCorr, 2);

    % Convert to LLA
    estDFT.lla = ecef2lla(estDFT.pos');

    % Log LLA Position    
    llaDF(i,:) = estDFT.lla;

end

hold on
geoplot(llaDF(:,1),llaDF(:,2),'*')
legend('Original','Dual Frequency IF','Location','northwest')

geobasemap satellite


