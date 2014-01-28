%localization_experiment   ITD-based sound source localization 
% 

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/01/22
%   ***********************************************************************

clear
close all
clc

% Add utilities
addpath Tools


%% ALGORITHM SETTINGS
% 
% 
% Depending on these settings, the localization models are calibrated such
% that the estimated ITD can be mapped to its corresponding sound source
% azimuth. If you keep these settings constant, the existing calibration
% files can be re-used, which reduces the runtime of the experiment.

% Reference sampling frequency in Hertz
fs = 48E3;

% Window size in seconds 
winSec = 20E-3;
  
% Lowest and highest center frequency in Hertz of the gammatone filterbank
fLowHz  = 100;
fHighHz = 12000;


%% ACOUSTIC SETTINGS
% 
% 
% Define number of competing speech sources
nSpeakers = 2;

% Minimum distance between competing sound sources in degree
minDistance = 5;

% *************************************************************************
% List of different rooms:
% *************************************************************************
% 'SURREY_A' 'SURREY_ROOM_A' 'SURREY_ROOM_B' 'SURREY_ROOM_C' 'SURREY_ROOM_D'
%  
%  RT60 = 0s  RT60 = 0.32s    RT60 = 0.47s    RT60 = 0.68s   RT60 = 0.89s
%  DRR  = inf DRR  = 6.09dB   DRR  = 5.31dB   DRR  = 8.82dB  DRR  = 6.12dB
% 

% rooms = {'SURREY_A' 'SURREY_ROOM_A' 'SURREY_ROOM_B' 'SURREY_ROOM_C' 'SURREY_ROOM_D'};
rooms = {'SURREY_A'};


%% EVALUATION SETTINGS
% 
% 
% Number of acoustic mixtures for each acoustic condition
nMixtures = 5; 

% Absolute error boundary in degree
thresDeg = 10;

% Visualize localization results 
bVisualize = true;

   
%% PERFORM CALIBRATION
% 
% 
% Frequency range of gammatone analysis
fRangeHz = [fLowHz fHighHz];

% Learn mapping between ITD and sound source azimuth
calibrate_ITD_Broadband(fs,winSec);
calibrate_ITD_Subband(fs,winSec,fRangeHz);


%% INITIALIZE PARAMETERS
%
%
% Reset internal states of random number generator. This allows to use
% different settings, while still obtaining the "same" random matrix with
% sound source positions.
rng(0);

% Azimuth range of sound source positions
azRange = (-90:5:90)';

% Audio path
pathAudio = [pwd,filesep,'Audio',filesep];

% Scan for audio files
audioFiles = listFiles(pathAudio,'*.wav');

% Number of different conditions
nRooms     = numel(rooms);
nSentences = numel(audioFiles);
nAzim      = numel(azRange);

% Allocate memory
[pc1,pc2,rmse1,rmse2] = deal(zeros(nMixtures,nRooms));

% Matrix with randomized sound source positions
azPos = azRange(round(1+(nAzim-1) * rand(nMixtures,nSpeakers)));

% Prepare visualization
if bVisualize
    % Figure handles
    fig1 = figure(99); fig2 = figure(100); arrange(fig1,fig2); 
end


%% MAIN LOOP OF THE LOCALIZATION EXPERIMENT
%
%
% Loop over number of acoustic mixtures
for ii = 1 : nMixtures
    
    % Randomly select "nSpeakers" sentences
    files = {audioFiles(round(1+(nSentences-1) * rand(nSpeakers,1))).name};
    
    % Enforce a "minDistance" spacing between all sources
    while any(diff(sort(azPos(ii,:),'ascend'),1,2) <= minDistance)
        % Revise random initialization
        azPos(ii,:) = azRange(round(1+(nAzim-1) * rand(1,nSpeakers)));
    end

    % Read audio signals
    audio = readAudio(files,fs);

    if bVisualize(1)
        % Encode original sound source positions for plotting purpose
        bVisualize = [bVisualize(1) azPos(ii,:)];
    end
    
    % Loop over number of different rooms
    for rr = 1 : nRooms
           
        % Spatialize audio signals using HRTF processing
        binaural = spatializeAudio(audio,fs,azPos(ii,:),rooms{rr});
        
        % Broadband estimation of the sound source azimuth
        azimEst1 = estimate_Azimuth_Broadband(binaural,fs,winSec,nSpeakers,bVisualize);
        
        % Subband estimation of the sound source azimuth
        azimEst2 = estimate_Azimuth_Subband(binaural,fs,winSec,fRangeHz,nSpeakers,bVisualize);
        
        % Evaluate localization performance
        [pc1(ii,rr),rmse1(ii,rr)] = evalPerformance(azPos(ii,:),azimEst1,thresDeg);
        [pc2(ii,rr),rmse2(ii,rr)] = evalPerformance(azPos(ii,:),azimEst2,thresDeg);
    end
    
    % Small delay to allow visual inspection of the figures
    if bVisualize(1); pause(1); end
    
    % Report progress
    fprintf('\nLocalization experiment: %.2f %%',100*ii/nMixtures);
end    
    
fprintf('\n\nBroadband approach\n');
fprintf('Percentage correct: %.2f %%\n%',mean(pc1,1));
fprintf('RMSE: %.2f %',nanmean(rmse1,1));
fprintf('\n\nSubband approach\n');
fprintf('Percentage correct: %.2f %%\n%',mean(pc2,1));
fprintf('RMSE: %.2f %',nanmean(rmse2,1));

