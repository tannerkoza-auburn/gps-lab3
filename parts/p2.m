clear; clc; close all

load('RCVR0.mat');
load('RCVRT.mat');



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
mean(errNormA)
std(errNormA)



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

mean(errNormB)
std(errNormB)

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
lambda = C/L1;
for i = 1:length(idx)
    
    % Novatel carrier in cycles.... trimble carrier in meters
    userPsr = RCVR0{idx(i,1),1}.L1.psr;
    userSats = RCVR0{idx(i,1),1}.L1.SVs;
    svPos = RCVR0{idx(i,1),1}.L1.svPos;
    userPhi = lambda*RCVR0{idx(i,1),1}.L1.carr; % convert carrier to m
    
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
legend('Novatel','Trimble')

errC = r_ab;
errNormC = vecnorm(errC);
time = nov.gpsTime(idx(:,1));
time = time - time(1);

mean(errNormC)
std(errNormC)

figure()
plot(time, errNormC, 'linewidth', 2)
title('Smoothed-Code Single-Difference Base Length Error')
ylabel('Error (m)')
xlabel('Time (s)')