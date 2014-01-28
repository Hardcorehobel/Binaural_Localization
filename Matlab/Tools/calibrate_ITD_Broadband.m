function calibrate_ITD_Broadband(fs,winSec,bCalibrate)
%calibrate_ITD_Broadband   Calculate ITD2azimuth mapping
%
%USAGE
%      calibrate_ITD_Broadband(fs,winSec)
%      calibrate_ITD_Broadband(fs,winSec,bCalibrate)
%
%INPUT PARAMETERS
%           fs : sampling frequency in Hertz
%       winSec : frame size in seconds of the cross-correlation analysis
%   bCalibrate : if true, enforce re-computation of the mapping function
% 
%OUTPUT PARAMETERS
%     The ITD2Azimuth mapping will be stored in the MAT file
%     ITD2Azimuth_Broadband.mat inside the \Tools directory.

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/01/22
%   v.0.2   2014/01/28 added flag to enforce re-calibration 
%   ***********************************************************************

% Initialize persistent variables
persistent PERfs PERwinSec


%% 1. CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 3 || isempty(bCalibrate); bCalibrate = false; end


%% 2. CALIBRATION SETTINGS
% 
% 
% Length of random noise sequence in seconds
lengthSec = 0.25;

% Use anechoic HRTF measurements
room = 'SURREY_A';

% Azimuth range of interest (real sound source positions)
azimRange = (-90:5:90);

% New azimuth range after interpolation
azimRangeInterp = (-90:1:90);

% Order of polynomial fit that is applied to the ITD to ensure a monotonic
% mapping
pOrder = 3;


%% 3. CALIBRATION STAGE
% 
% 
% Check if we can re-use the calibration file from the last function call
if isequal(fs,PERfs) && isequal(winSec,PERwinSec) && ~bCalibrate
    bRecalibrate = false;
else
    bRecalibrate = true;
end

% Perform calibration
if bRecalibrate
    % Store persistent variables
    PERfs = fs; PERwinSec = winSec;
    
    % Number of different sound source positions
    nAzim = numel(azimRange);
    
    % Create white noise
    noise = randn(round(lengthSec*fs),1);
    
    % Allocate memory
    itd2Azim = zeros(nAzim,1);
    
    
    % MAIN LOOP
    %
    %
    % Loop over number of different sound source directions
    for ii = 1 : nAzim
        
        % Spatialize audio signal using HRTF processing
        binaural = spatializeAudio(noise,fs,azimRange(ii),room);

        % Estimate ITD
        [itdEst,lags] = estimate_ITD_Broadband(binaural,fs,winSec);
        
        % Store azimuth-dependent ITD
        itd2Azim(ii) = itdEst;
        
        % Report progress
        fprintf('\nBroadband-based ITD2Azimuth calibration: %.2f %%',100*ii/nAzim);
    end
    
    
    % Interpolation
    %
    %
    % Interpolate to 'rangeAzInterp'
    itd2AzimInterp = interp1(azimRange,itd2Azim,azimRangeInterp);
    
    % Ensure that mapping is monotonic by using a polynomial fit
    itd2AzimPoly = polyval(polyfit(azimRangeInterp,itd2AzimInterp,pOrder),azimRangeInterp);
        
    
    % Save data
    %
    %
    mapping.fs           = fs;
    mapping.azim         = azimRangeInterp;
    mapping.itd          = lags/fs;
    mapping.itd2azimRaw  = itd2Azim;
    mapping.itd2azim     = itd2AzimPoly;
    mapping.polyOrder    = pOrder;
    mapping.itdMax       = max(itd2AzimPoly);
    mapping.itdMin       = min(itd2AzimPoly);
    
    % Store ITD2Azimuth template
    save([pwd,filesep,'Tools',filesep,'ITD2Azimuth_Broadband.mat'],'mapping');
end