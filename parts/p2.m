clear; clc; close all

load('RCVR0.mat');
load('RCVRT.mat');

% truth given from surveyed position of trimble antenna
basePosTruth = [423203.359; -5361678.541; 3417280.681];



%% ---- part a

% instantiate classes
novClass = gnssReceiver(0.15); % uses RCVR0 data
trimClass = gnssReceiver(0.15); % uses RCVRT data

% --- solve Novatel
for i = 1:length(RCVR0)
    psr = RCVR0{i,1}.L1.psr;
    dopp = RCVR0{i,1}.L1.dopp;
    svPos = RCVR0{i,1}.L1.svPos;
    svVel = RCVR0{i,1}.L1.svVel;
    svClockCorr = RCVR0{i,1}.L1.clkCorr;
    carrFreq = 1;
    novatel{i,1} = novClass.pv3D(psr, dopp, svPos, svVel, svClockCorr, carrFreq);
    novatel{i,1}.gpsTime = RCVR0{i,1}.L1.gpsTime;
    
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
geoplot(trim.poslla(:,1), trim.poslla(:,2), '*')
geobasemap satellite
title('Receiver Position Solutions')
legend('Novatel','Trimble')


% --- sync solutions in time
% novatel starts first --- find starting novatel which may have a match
[M,strt] = min(abs(nov.gpsTime - trim.gpsTime(1)));

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
errA = nov.pos(:, idx(:,1)) - trim.pos(:, idx(:,2));
errNormA = vecnorm(nov.pos(:, idx(:,1)) - trim.pos(:, idx(:,2)));


% err
figure()
subplot(3,1,1)
plot(time, errA(1,:), 'linewidth', 2);
subplot(3,1,2)
plot(time, errA(2,:), 'linewidth', 2);
subplot(3,1,3)
plot(time, errA(3,:), 'linewidth', 2);

figure()
plot(time,errNormA, 'linewidth', 2);
title('Trimble and Novatel Base Length Error')
ylabel('Error (m)')
xlabel('Time (s)')
meanA = mean(errNormA)
stdA = std(errNormA)



%% ---- part B

% --- used sync data steps to compute DGPS solution
for i = 1:length(idx)
    userPsr = RCVR0{idx(i,1),1}.L1.psr;
    userSats = RCVR0{idx(i,1),1}.L1.SVs;
    svPos = RCVR0{idx(i,1),1}.L1.svPos;
    
    basePsr = RCVRT{idx(i,2),1}.L1.psr;
    basePos = trim.pos(:,idx(i,2));
    baseSats = RCVRT{idx(i,2),1}.L1.SVs;
    
    % find matching satellites at this instance in time
    [C, iUser, iBase] = intersect(userSats, baseSats);
   
    % --- trim to only use matching pseudoranges
    userPsr = userPsr(iUser);
    svPos = svPos(iUser,:);
    basePsr = basePsr(iBase);
    
    [out] = novClass.sdp3D(userPsr,basePsr,svPos, basePos);
    
    r_ab(:,i) = basePos - out.pos; % find relative position    
end 
errB = r_ab;
errNormB = vecnorm(errB);
time = nov.gpsTime(idx(:,1));
time = time - time(1);

meanB = mean(errNormB)
stdB = std(errNormB)

figure()
plot(time, errNormB, 'linewidth', 2)
title('Single-Difference Base Length Error')
ylabel('Error (m)')
xlabel('Time (s)')



%% --- part C

% --- differential smoothing
M = length(idx);
C = physconst('LightSpeed');
L1 = 1575.42 * 10^6; % freq of L1
L2 = 1227.6 * 10^6; % freq of L2
lambdaL1 = C/L1;
lambdaL2 = C/L2;
for i = 1:length(idx)
    
    % Novatel carrier in cycles.... trimble carrier in meters
    userPsr = RCVR0{idx(i,1),1}.L1.psr;
    userSats = RCVR0{idx(i,1),1}.L1.SVs;
    svPos = RCVR0{idx(i,1),1}.L1.svPos;
    userPhi = lambdaL1*RCVR0{idx(i,1),1}.L1.carr; % convert carrier to m
    
    basePsr = RCVRT{idx(i,2),1}.L1.psr;
    basePos = trim.pos(:,idx(i,2));
    baseSats = RCVRT{idx(i,2),1}.L1.SVs;
    basePhi = RCVRT{idx(i,2),1}.L1.carr;
    
    % find matching satellites at this instance in time
    [C, iUser, iBase] = intersect(userSats, baseSats);
    
    
    % --- trim to only use matching pseudoranges
    userPsr = userPsr(iUser);
    userPhi = userPhi(iUser);
    svPos = svPos(iUser,:);
    basePsr = basePsr(iBase);
    basePhi = basePhi(iBase);
    
    
    % smooth each pseudorange 
    for j = 1:length(basePsr) 
        
        % i = time j = psr
        delRho(j,i) = userPsr(j) - basePsr(j);
        delPhi(j,i) = userPhi(j) - basePhi(j);
        
        if i == 1
            delRhoBar(j,i) = delRho(j,i);
        else
            delRhoBar(j,i) = (1/M)*delRho(j,i) + ((M-1)/M)*(delRhoBar(j,i-1) + delPhi(j,i) - delPhi(j,i-1));
        end
        
    end
    
    [out] = novClass.sdp3D(delRhoBar(:,i), 0*delRhoBar(:,i),svPos, basePos);
    
    smooth_DGPS(:,i) = out.pos;
    
    r_ab(:,i) = basePos - out.pos; % find relative position    
       
end 


% --- plot solutions
smooth_lla = ecef2lla(smooth_DGPS');
figure()
geoplot(smooth_lla(:,1), smooth_lla(:,2), '*')
geobasemap satellite
title('Smoothed DGPS Novatel Solution')

errC = r_ab;
errNormC = vecnorm(errC);
time = nov.gpsTime(idx(:,1));
time = time - time(1);

meanC = mean(errNormC)
stdC = std(errNormC)

figure()
plot(time, errNormC, 'linewidth', 2)
title('Smoothed-Code Single-Difference Base Length Error')
ylabel('Error (m)')
xlabel('Time (s)')




%% --- part D RTK

%--- solve for the integer ambiguities


% unpack psr, carrL1 and carrL2 for each sat... only performing once

SVintsL1 = NaN(32,length(idx));
SVintsL2 = NaN(32,length(idx));
delPsr_CarrL1 = NaN(32, length(idx));

for i = 1:length(idx)
    % Novatel carrier in cycles.... trimble carrier in meters
    userPsr = RCVR0{idx(i,1),1}.L1.psr;
    userCarrL1 = lambdaL1*RCVR0{idx(i,1),1}.L1.carr;
    userSatsL1 = RCVR0{idx(i,1),1}.L1.SVs;
    userCarrL2 = lambdaL2*RCVR0{idx(i,1),1}.L2.carr;
    userSatsL2 = RCVR0{idx(i,1),1}.L2.SVs;
    svPos = RCVR0{idx(i,1),1}.L1.svPos;
     %userPhi = lambdaL1*RCVR0{idx(i,1),1}.L1.carr; % convert carrier to m
    
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
    basePos = basePosTruth; % use true base position?
    
    
      
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
    err_d_sd(i) = norm(out.rpv); % single difference storage
    err_d_dd(i) = norm(outDD.rpv); % double difference storage
    
end
meanD_sd = mean(err_d_sd)
stdD_sd = std(err_d_sd)
meanD_sd = mean(err_d_dd)
stdD_sd = std(err_d_dd)



%--- plot SV int over time
figure()
subplot(2,1,1)
plot(time, SVintsL1');
ylim([-200 200])
title('L1 Single-Difference Integer Ambiguity Estimates')
ylabel('Integer')
xlabel('Time (s)')
subplot(2,1,2)
plot(time, SVintsL2');
ylim([-200 200])
title('L2 Single-Difference Integer Ambiguity Estimates')
ylabel('Integer')
xlabel('Time (s)')

figure()
plot(time, err_d_sd, 'lineWidth', 2);
title('L1 Carrier-Based Base Length Error (Widelane single diff)')
ylabel('Error (m)')
xlabel('Time (s)')

figure()
plot(time, err_d_dd, 'lineWidth', 2);
title('L1 Carrier-Based Base Length Error (Widelane double diff)')
ylabel('Error (m)')
xlabel('Time (s)')

% ------ perform solution again using averaged integers over whole sample window
meanInts = mean(delPsr_CarrL1, 2, 'omitnan');
for i = 1:length(idx)
    % Novatel carrier in cycles.... trimble carrier in meters
    userPsr = RCVR0{idx(i,1),1}.L1.psr;
    userCarrL1 = lambdaL1*RCVR0{idx(i,1),1}.L1.carr;
    userSatsL1 = RCVR0{idx(i,1),1}.L1.SVs;
    userCarrL2 = lambdaL2*RCVR0{idx(i,1),1}.L2.carr;
    userSatsL2 = RCVR0{idx(i,1),1}.L2.SVs;
    svPos = RCVR0{idx(i,1),1}.L1.svPos;
     %userPhi = lambdaL1*RCVR0{idx(i,1),1}.L1.carr; % convert carrier to m
    
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
    basePos = basePosTruth; % use true base position?
    
    
      
    % use solved-for integers in the single-difference algorithm:
    userPsr = userCarrL1 - baseCarrL1;
    basePsr = -lambdaL1*round(meanInts(matchingSats)); % corresponding integers to available sats
    out = novClass.sdp3D(userPsr, basePsr, svPos, basePos);
    
   
    err_d(i) = norm(out.rpv);
    
end
meanD = mean(err_d)
stdD = std(err_d)


figure()
plot(time, err_d, 'lineWidth', 2);
title('L1 Carrier-Based Base Length Error (Averaging Window)')
ylabel('Error (m)')
xlabel('Time (s)')


