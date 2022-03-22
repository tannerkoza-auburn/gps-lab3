%% Fundamentals of GPS - Lab Data Parser

clear
clc

%% Restructure Data

% Directory Declaration
indata = dir(fullfile('data/receiver_data','*.mat'));
outdir = 'data/parsed_data';

% Define Number of Files
numFiles = length(indata);

if numFiles == 0
    error('No files have been read. Run this file from the root directory.')
end

% Define Constants
C = physconst('LightSpeed');

for i = 1:numFiles
    thisfile = indata(i).name;
    outfile = strrep(thisfile,'_data','');
    outstruct = strrep(outfile,'.mat','');
    destfile = fullfile(outdir,outfile);

    tstruct = load(thisfile);
    structname = fieldnames(tstruct);
    sstruct = structname{1};

    % Preallocation
    outdata = cell(length(tstruct.(sstruct).gpsTime(:,1)),1);

    % Parsing & Data Insertion
    for j = 1:length(tstruct.(sstruct).gpsTime(:,1))
        
        % Find SVs in View
        [~, satsInView] = find(~isnan(tstruct.(sstruct).obs.L1.psr(j,:)));

        % Extract PSR, Dopp, and Time
        outdata{j}.gpsTime = tstruct.(sstruct).gpsTime(j,1);
        outdata{j}.psr(:,1) = tstruct.(sstruct).obs.L1.psr(j,satsInView);
        outdata{j}.dopp(:,1) = tstruct.(sstruct).obs.L1.dopp(j,satsInView);

        % Calculate SV Position and Velocity
        transitTime = outdata{j}.psr/C;
        transmitTime = tstruct.(sstruct).gpsTime(j) - transitTime;

        for k = 1:length(transmitTime)
            [outdata{j}.pos(k,:), outdata{j}.vel(k,:), outdata{j}.clkCorr(k,1)] = calc_sv_pos(tstruct.(sstruct).ephem(satsInView(k),:), transmitTime(k), transitTime(k));
        end 

    end

    newstruct.(outstruct) = outdata;
    save(destfile, '-struct', 'newstruct')

    clearvars outdata newstruct

end
