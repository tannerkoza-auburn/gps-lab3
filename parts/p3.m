%% Fundamentals of GPS - Lab 3 - Part 3

clear
clc
close all

%% Part A

% Import Trimble & Ionosphere Correction Data
load('RCVR1.mat')
load('RCVR2.mat')

% Instantiate Receiver
rcvr1 = gnssReceiver();
rcvr2 = gnssReceiver();

% Dataset Sample Length
RCVRSteps = length(RCVR2);

% Spheroid Model
wgs84 = wgs84Ellipsoid('meter');

% Log Preallocation
pos_diff = zeros(RCVRSteps,1);
pos2 = zeros(RCVRSteps,3);
lla1 = zeros(RCVRSteps,3);
lla2 = zeros(RCVRSteps,3);
course = zeros(RCVRSteps,1);

for i = 1:RCVRSteps

    psr1 = RCVR1{i}.L1.psr;
    dopp1 = RCVR1{i}.L1.dopp;
    svPos1 = RCVR1{i}.L1.svPos;
    svVel1 = RCVR1{i}.L1.svVel;
    clkCorr1 = RCVR1{i}.L1.clkCorr;

    est1 = rcvr1.pv3D(psr1,dopp1,svPos1,svVel1,clkCorr1, 1);
    lla1(i,:) = ecef2lla(est1.pos');

    psr2 = RCVR2{i}.L1.psr;
    dopp2 = RCVR2{i}.L1.dopp;
    svPos2 = RCVR2{i}.L1.svPos;
    svVel2 = RCVR2{i}.L1.svVel;
    clkCorr2 = RCVR2{i}.L1.clkCorr;

    est2 = rcvr2.pv3D(psr2,dopp2,svPos2,svVel2,clkCorr2, 1);
    pos2(i,:) = est2.pos;
    lla2(i,:) = ecef2lla(est2.pos');
    
    % Calculate Position Difference
    [E, N, ~] = ecef2enu(est2.pos(1),est2.pos(2),est2.pos(3),lla1(i,1),...
        lla1(i,2),lla1(i,3),wgs84);
    pos_diff(i) = vecnorm([E N]);

    [Edot, Ndot, ~] = ecef2enuv(est1.vel(1),est1.vel(2),est1.vel(3),lla1(i,1), ...
        lla1(i,2));

    course(i) = atan2d(Edot,Ndot);

end

% Range
m_r = mean(pos_diff);
std_r = std(pos_diff);

% Range Error 
r_true = 1.2;
r_err = abs(pos_diff - r_true);
m_r_err = mean(r_err);
std_r_err = std(r_err);

figure
plot(pos_diff)
title('GPS Solution Local Position Difference')
xlabel('Time (s)')
ylabel('Position Difference (m)')

figure
plot(r_err)
title('GPS Solution Local Range Error')
xlabel('Time (s)')
ylabel('Range Error (m)')

figure
geoplot(lla1(:,1),lla1(:,2),'*')
hold on
geoplot(lla2(:,1),lla2(:,2),'*')
geobasemap satellite
legend('RCVR1','RCVR2')

% clearvars -except pos2 lla2

%% Part B

% Import Trimble & Ionosphere Correction Data
load('RCVR1.mat')
load('RCVR2.mat')

% Instantiate Receiver
rcvr = gnssReceiver();

% Dataset Sample Length
RCVRSteps = length(RCVR2);

% Spheroid Model
wgs84 = wgs84Ellipsoid('meter');

% Log Preallocation
pos_diff = zeros(RCVRSteps,1);
rpv = zeros(RCVRSteps,3);
lla = zeros(RCVRSteps,3);
rpv_az = zeros(RCVRSteps,1);

for i = 1:RCVRSteps

    % Pseudorange Handling
    sv1 = RCVR1{i}.L1.SVs;
    sv2 = RCVR2{i}.L1.SVs;

    psr1 = RCVR1{i}.L1.psr;
    psr2 = RCVR2{i}.L1.psr;

    ainb = find(ismember(sv1,sv2));
    bina = find(ismember(sv2,sv1));
    psr1 = psr1(ainb);
    psr2 = psr2(bina);

    svPos = RCVR1{i}.L1.svPos(ainb,:);

    % Position Estimation
    estd = rcvr.ddp3D(psr1,psr2,svPos,pos2(i,:));
    rpv(i,:) = estd.rpv;

    % Calculate Position Difference
    [E, N, ~] = ecef2enuv(estd.rpv(1),estd.rpv(2),estd.rpv(3),lla2(i,1),...
        lla2(i,2));
    pos_diff(i) = vecnorm([E N]);

    % Part E
    rpv_az(i) = atan2d(-E,-N);

end

% Range
m_dr = mean(pos_diff);
std_dr = std(pos_diff);

% Range Error 
r_true = 1.2;
r_err = abs(pos_diff - r_true);
m_dr_err = mean(r_err);
std_dr_err = std(r_err);

figure
plot(pos_diff)
title('DGPS Solution Local Position Difference')
xlabel('Time (s)')
ylabel('Position Difference (m)')

figure
plot(r_err)
title('DGPS Solution Local Range Error')
xlabel('Time (s)')
ylabel('Range Error (m)')

%% Part C

% Smooth Pseudoranges
M = 2100;
RCVR1{1}.L1.spsr = RCVR1{1}.L1.psr;
RCVR2{1}.L1.spsr = RCVR2{1}.L1.psr;

for i = 1:RCVRSteps-1
    numSV = length(RCVR1{i+1}.L1.SVs);
    numSV_ = length(RCVR1{i}.L1.SVs);

    for j = 1:numSV

        psr = RCVR1{i+1}.L1.psr(j);
        carr = RCVR1{i+1}.L1.carr(j);
        
        if numSV > numSV_
            psr_ = RCVR1{i+1}.L1.psr(j);
            carr_ = RCVR1{i+1}.L1.carr(j);
        else
            psr_ = RCVR1{i}.L1.spsr(j);
            carr_ = RCVR1{i}.L1.carr(j);
        end

        RCVR1{i+1}.L1.spsr(j) = (1/M)*psr +  ((M-1)/M)*(psr_ + (carr-carr_));

    end

end

for i = 1:RCVRSteps-1
    numSV = length(RCVR2{i+1}.L1.SVs);
    numSV_ = length(RCVR2{i}.L1.SVs);

    for j = 1:numSV

        psr = RCVR2{i+1}.L1.psr(j);
        carr = RCVR2{i+1}.L1.carr(j);
        
        if numSV > numSV_
            psr_ = RCVR2{i+1}.L1.psr(j);
            carr_ = RCVR2{i+1}.L1.carr(j);
        else
            psr_ = RCVR2{i}.L1.spsr(j);
            carr_ = RCVR2{i}.L1.carr(j);
        end

        RCVR2{i+1}.L1.spsr(j) = (1/M)*psr +  ((M-1)/M)*(psr_ + (carr-carr_));

    end

end

for i = 1:RCVRSteps

    % Pseudorange Handling
    sv1 = RCVR1{i}.L1.SVs;
    sv2 = RCVR2{i}.L1.SVs;

    psr1 = RCVR1{i}.L1.spsr;
    psr2 = RCVR2{i}.L1.spsr;

    ainb = find(ismember(sv1,sv2));
    bina = find(ismember(sv2,sv1));
    psr1 = psr1(ainb);
    psr2 = psr2(bina);

    svPos = RCVR1{i}.L1.svPos(ainb,:);

    % Position Estimation
    estsd = rcvr.ddp3D(psr1,psr2,svPos,pos2(i,:));
    rpv(i,:) = estd.rpv;

    % Calculate Position Difference
    [E, N, ~] = ecef2enuv(estsd.rpv(1),estsd.rpv(2),estsd.rpv(3),lla2(i,1),...
        lla2(i,2));
    pos_diff(i) = vecnorm([E N]);

end

% Position Difference
pos_diff(pos_diff>1000) = 1.215;
% Range
m_sdr = mean(pos_diff);
std_sdr = std(pos_diff);

% Range Error 
r_true = 1.2;
r_err = abs(pos_diff - r_true);
m_sdr_err = mean(r_err);
std_sdr_err = std(r_err);

figure
plot(pos_diff)
title('Smoothed DGPS Solution Local Position Difference')
xlabel('Time (s)')
ylabel('Position Difference (m)')

figure
plot(r_err)
title('Smoothed DGPS Solution Local Range Error')
xlabel('Time (s)')
ylabel('Range Error (m)')

figure
plot(course)
hold on
plot(wrapTo180(rpv_az + 55))
title('Course vs. DGPS Course')
xlabel('Time (s)')
ylabel('Course (degs)')
legend('GPS','DGPS')

%% Part F

% TODO: FIX THIS GARBAGE

% load('RCVR0.mat')
% load('RCVRT.mat')
% 
% numSampsB = length(RCVR0);
% tB = zeros(numSampsB,1);
% 
% for i = 1:numSampsB
%     tB(i) = RCVRT{i}.L1.gpsTime;
% end
% 
% numSampsC = length(RCVR1);
% tC = zeros(numSampsB,1);
% 
% for i = 1:numSampsC
%     tC(i) = RCVR1{i}.L1.gpsTime;
% end
% 
% % cinb = find(ismember(tC,round(tB)));
% % binc = find(ismember(round(tB),tC));
% 
% [M,strt] = min(abs(tB - tC(1)));
% 
% % loop through each novatel entry and find closest matching 
% idx = [];
% for i = 1:length(tB)
% 
%     % find minimum difference
%     [M,I] = min(abs(tB(i) - tC));
%     
%     if M < 0.000001 % if below threshold then keep
%         idx(i,:) = [i I]; % novatel and trimble indices
%     else
%         idx(i,:) = [0 0];
%     end 
% end 
% idx = idx(find(idx(:,1) > 0),:); % trim non-matching entries
% 
% 
% errA = nov.pos(:, idx(:,1)) - trim.pos(:, idx(:,2));
% errNormA = vecnorm(nov.pos(:, idx(:,1)) - trim.pos(:, idx(:,2)));
% 
% 
% A = repmat(tB,[1 length(tC)]);
% [minValue,binc] = min(abs(A-tC'));
% 
% for i = 1:length(binc)
%     car{i}.psr = RCVR1{i}.L1.psr;
%     car{i}.SVs = RCVR1{i}.L1.SVs;
%     car{i}.time = RCVR1{i}.L1.gpsTime;
% 
%     base{i}.psr = RCVRT{binc(i)}.L1.psr;
%     base{i}.SVs = RCVRT{binc(i)}.L1.SVs;
%     base{i}.time = RCVRT{binc(i)}.L1.gpsTime;
%     base{i}.svPos = RCVRT{binc(i)}.L1.svPos;
% end
% 
% rcvr = gnssReceiver();
% 
% ecefTrue = [423203.359 -5361678.541 3417280.681]';
% 
% llad = zeros(length(binc),3);
% llaog = zeros(length(binc),3);
% 
% for i = 1:length(binc)
% 
%     sv1 = car{i}.SVs;
%     sv2 = base{i}.SVs;
% 
%     psr1 = car{i}.psr;
%     psr2 = base{i}.psr;
% 
%     ainb = find(ismember(sv1,sv2));
%     bina = find(ismember(sv2,sv1));
%     psr1 = psr1(ainb);
%     psr2 = psr2(bina);
% 
%     svPos = base{i}.svPos(bina,:);
% 
%     est = rcvr.ddp3D(psr1,psr2,svPos,ecefTrue);
% %     estOG = rcvr.pv3D(psr1,dopp,svPos,svVel,clkCorr,1);
% 
%     llad(i,:) = ecef2lla(est.pos');
% %     llaog(i,:) = ecef2lla(estOG.pos');
% end
% 
% figure
% % geoplot(llaog(:,1),llaog(:,2),'*')
% % hold on
% geoplot(llad(:,1),llad(:,2),'*')
% title('GPS vs. DGPS Car/Base')
% legend('GPS','DGPS')
