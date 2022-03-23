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

    % Frequency Determination
    freq = fieldnames(tstruct.(sstruct).obs);
    numFreq = length(freq);

    % Preallocation
    outdata = cell(length(tstruct.(sstruct).gpsTime(:,1)),1);

    % Parsing & Data Insertion
    for j = 1:length(tstruct.(sstruct).gpsTime(:,1))
        
        for f = 1:numFreq

            % Find SVs in View
            [~, satsInView] = find(~isnan(tstruct.(sstruct).obs.(freq{f}).psr(j,:)));
    
            % Extract SVs in View, PSR, Dopp, and Time
            outdata{j}.(freq{f}).SVs = satsInView;
            outdata{j}.(freq{f}).gpsTime = tstruct.(sstruct).gpsTime(j,1);
            outdata{j}.(freq{f}).psr(:,1) = tstruct.(sstruct).obs.(freq{f}).psr(j,satsInView);
            outdata{j}.(freq{f}).dopp(:,1) = tstruct.(sstruct).obs.(freq{f}).dopp(j,satsInView);
            outdata{j}.(freq{f}).carr(:,1) = tstruct.(sstruct).obs.(freq{f}).carr(j,satsInView);
    
            % Calculate SV Position and Velocity
            transitTime = outdata{j}.(freq{f}).psr/C;
            transmitTime = tstruct.(sstruct).gpsTime(j) - transitTime;

            for k = 1:length(transmitTime)
                [outdata{j}.(freq{f}).svPos(k,:), outdata{j}.(freq{f}).svVel(k,:), outdata{j}.(freq{f}).clkCorr(k,1)] = calc_sv_pos(tstruct.(sstruct).ephem(satsInView(k),:), transmitTime(k), transitTime(k));
            end 

        end
    


        clearvars satsInView

    end

    newstruct.(outstruct) = outdata;
    save(destfile, '-struct', 'newstruct')

    clearvars outdata newstruct

end
