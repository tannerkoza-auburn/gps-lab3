%{
    DGPS with dynamic base length. Trimble on roof of Woltosz and Novatel
    on roof of car. Will attempt Lambda-method RTK and compare to
    traditional single-difference DGPS solution (code-only).
%}

clear; clc; close all;

%% --- Import Data
load('RCVRT.mat'); % trimble (static on woltosz roof)
load('RCVR1.mat'); % novatel on car

% truth given from surveyed position of trimble antenna
basePosTruth = [423203.359; -5361678.541; 3417280.681];


%% ---- Preliminary baseline solution (NO DGPS)

% instantiate classes
novClass = gnssReceiver(0.15); % uses RCVR1 data
trimClass = gnssReceiver(0.15); % uses RCVRT data

% --- solve Novatel
for i = 1:length(RCVR1)
    psr = RCVR1{i,1}.L1.psr;
    dopp = RCVR1{i,1}.L1.dopp;
    svPos = RCVR1{i,1}.L1.svPos;
    svVel = RCVR1{i,1}.L1.svVel;
    svClockCorr = RCVR1{i,1}.L1.clkCorr;
    carrFreq = 1;
    novatel{i,1} = novClass.pv3D(psr, dopp, svPos, svVel, svClockCorr, carrFreq);
    novatel{i,1}.gpsTime = RCVR1{i,1}.L1.gpsTime;
    
    % extraction for plotting
    nov.pos(:,i) = novatel{i,1}.pos;
    nov.gpsTime(i) = novatel{i,1}.gpsTime;
end 

% --- solve Trimble
for i = 1:length(RCVRT)    
    psr = RCVRT{i,1}.L1.psr;
    dopp = RCVRT{i,1}.L1.dopp;
    svPos = RCVRT{i,1}.L1.svPos;
    svVel = RCVRT{i,1}.L1.svVel;
    svClockCorr = RCVRT{i,1}.L1.clkCorr;
    carrFreq = 1;
    trimble{i,1} = trimClass.pv3D(psr, dopp, svPos, svVel, svClockCorr, carrFreq);    
    trimble{i,1}.gpsTime = RCVRT{i,1}.L1.gpsTime;
    
    % extraction for plotting 
    trim.pos(:,i) = trimble{i,1}.pos;
    trim.gpsTime(:,i) = trimble{i,1}.gpsTime;    
end 


% --- plot solutions
nov.poslla = ecef2lla(nov.pos');
trim.poslla = ecef2lla(trim.pos');
geoplot(nov.poslla(:,1), nov.poslla(:,2), '*')
hold on
geoplot(trim.poslla(1,1), trim.poslla(1,2), '*', 'MarkerSize', 15)
geobasemap satellite
title('Receiver Position Solutions')
legend('Novatel','Trimble')


%% --- Sync data sets in time
%{
    idx contains corresponding indices for novatel and trimble as column
    vectors, respectively. eg: idx(1,1) and idx(1,2) are the indices in
    trimble and novatel vectors which correspond to each other in time
%}

% loop through each novatel entry and find closest matching 
idx = [];
for i = 1:length(nov.gpsTime)

    % find minimum difference
    [M,I] = min(abs(nov.gpsTime(i) - trim.gpsTime));
    
    if M < 0.000001 % if below threshold then keep
        idx(i,:) = [i I]; % novatel and trimble indices
    else
        idx(i,:) = [0 0];
    end 
end 
idx = idx(find(idx(:,1) > 0),:); % trim non-matching entries 
time = nov.gpsTime(idx(:,1));
time = time - time(1); % reset time



%% --- DGPS solutions


C = physconst('LightSpeed');
L1 = 1575.42 * 10^6; % freq of L1
L2 = 1227.6 * 10^6; % freq of L2
lambdaL1 = C/L1;
lambdaL2 = C/L2;

%--- solve for the integer ambiguities
% unpack psr, carrL1 and carrL2 for each sat... 

for i = 1:length(idx)
    % Novatel carrier in cycles.... trimble carrier in meters
    userPsr = RCVR1{idx(i,1),1}.L1.psr;
    userCarrL1 = RCVR1{idx(i,1),1}.L1.carr;
    userSatsL1 = RCVR1{idx(i,1),1}.L1.SVs;
    userCarrL2 = RCVR1{idx(i,1),1}.L2.carr;
    userSatsL2 = RCVR1{idx(i,1),1}.L2.SVs;
    svPos = RCVR1{idx(i,1),1}.L1.svPos;
     %userPhi = lambdaL1*RCVR1{idx(i,1),1}.L1.carr; % convert carrier to m
    
    basePos = trim.pos(:,idx(i,2)); % position used for DGPS base station
    basePsr = RCVRT{idx(i,2),1}.L1.psr;
    baseCarrL1 = RCVRT{idx(i,2),1}.L1.carr;
    baseSatsL1 = RCVRT{idx(i,2),1}.L1.SVs;
    baseCarrL2 = RCVRT{idx(i,2),1}.L2.carr;
    baseSatsL2 = RCVRT{idx(i,2),1}.L2.SVs;
    
    % find matching satellites for btwn base and user for L1 and L2
    [C, ~, ~] = intersect(userSatsL1, userSatsL2); % user L1 and L2
    [C, ~, ~] = intersect(baseSatsL1, C);  % base L1 and resulting
    % find matching over all
    [matchingSats, ~, ~] = intersect(baseSatsL2, C); % resulting and base L2
    
    % now extract all matching satellite indices
    [~, iUserL1, ~] = intersect(userSatsL1, matchingSats);
    [~, iUserL2, ~] = intersect(userSatsL2, matchingSats);
    [~, iBaseL1, ~] = intersect(baseSatsL1, matchingSats);
    [~, iBaseL2, ~] = intersect(baseSatsL2, matchingSats);
    
    
    % --- trim to only use matching pseudoranges, carr, satPos
    userPsr = userPsr(iUserL1);
    userCarrL1 = userCarrL1(iUserL1);
    userCarrL2 = userCarrL2(iUserL2);
    svPos = svPos(iUserL1,:);
    basePsr = basePsr(iBaseL1);
    baseCarrL1 = baseCarrL1(iBaseL1);
    baseCarrL2 = baseCarrL2(iBaseL2);
    basePos = basePosTruth;
    
    % push through single-difference carrier based solution (solves for integer ambiguities)
    out = novClass.sdCarr3D(userPsr, userCarrL1, userCarrL2, basePsr, baseCarrL1, baseCarrL2, svPos, basePos);
    outDD = novClass.ddCarr3D(userPsr, userCarrL1, userCarrL2, basePsr, baseCarrL1, baseCarrL2, svPos, basePos);
    
     % save integer ambiguities with associated SV #
    for j = 1:length(matchingSats)
        intIdx = matchingSats(j); % which sat?
        SVintsL1(intIdx,i) = out.N_lambda(2*j-1);
        SVintsL2(intIdx,i) = out.N_lambda(2*j);
        
        % storing deltas for estimating integers w/ averaging window
        delPsr_CarrL1(intIdx, i) = ((userPsr(j) - basePsr(j)) - (userCarrL1(j) - baseCarrL1(j)))/lambdaL1;
    end     
    
    
    % store position solutions
    nov.posDGPS(:,i) = out.pos; % single diff carrier
    nov.posDGPS_DD(:,i) = outDD.pos; % double-diff carrier
end

% --- plot solutions on top of standrd GPS
nov.posllaDGPS = ecef2lla(nov.posDGPS');
nov.posllaDGPS_DD = ecef2lla(nov.posDGPS_DD');
geoplot(nov.posllaDGPS(:,1), nov.posllaDGPS(:,2), '*')
geoplot(nov.posllaDGPS_DD(:,1), nov.posllaDGPS_DD(:,2), '*')
geobasemap satellite
title('Receiver Position Solutions')
legend('GPS','Base Station', 'DGPS_{\Delta \Phi}', 'DGPS_{\Delta \nabla \Phi}')


figure()
subplot(2,1,1)
plot(SVintsL1');
ylim([-200 200])
title('L1 Single-Difference Integer Ambiguity Estimates')
ylabel('Integer')
xlabel('Time (s)')
subplot(2,1,2)
plot(SVintsL2');
ylim([-200 200])
title('L2 Single-Difference Integer Ambiguity Estimates')
ylabel('Integer')
xlabel('Time (s)')

