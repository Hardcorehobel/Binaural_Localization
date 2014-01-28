function [itd,lags] = estimate_ITD_Subband(binaural,fs,winSec,fRange,N)
%estimate_ITD_Subband   Subband ITD estimation.
%
% The interaural time difference (ITD) is estimated from binaural signals.
% A peripheral auditory processing stage is used to decomposes the input 
% signals into individual frequency channels using a gammatone filterbank
% (gammaFB.m). The center frequencies are equally spaced on the equivalent
% rectangular bandwidth (ERB)-rate scale between fRange(1) and fRange(2).
% Then, the cross-correlation function (CCF) is computed for each subband
% over short time frames (winSec). The resulting 3-dimensional
% cross-correlation function, which is a function of the number of lags,
% subbands and frames (nLags x nSubbands x nFrames), is integrated
% across time frames and the most prominent peaks are assumed to reflect
% the estimated ITDs in each subband for N sources.          
% 
%USAGE
%    [itd,lags] = estimate_ITD_Subband(binaural,fs)
%    [itd,lags] = estimate_ITD_Subband(binaural,fs,winSec,fRange,N)
%
%INPUT PARAMETERS
%   binaural : binaural input signal [nSamples x 2]
%         fs : sampling frequency in Hertz
%     winSec : frame size in seconds of the cross-correlation analysis
%              (default, winSec = 20E-3)
%     fRange : [fLowHz fHighHz] where fLowHz and fHighHz define the lowest
%              and the highest center frequency of the gammatone filterbank
%              (default, fRange = [100 12000])
%          N : number of sources that should be localized (default, N = 1) 
%
%OUTPUT PARAMETERS
%        itd : subband-based ITD estimates of all N sources [nSubbands x N]
%       lags : time lags over which the CCF is computed

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/01/26
%   ***********************************************************************


%% 1. CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameters
if nargin < 3 || isempty(winSec);   winSec   = 20E-3;       end
if nargin < 4 || isempty(fRange);   fRange   = [100 12000]; end
if nargin < 5 || isempty(N);        N        = 1;           end


%% 2. INITIALIZE PARAMETERS
% 
% 
% Maximum time delay that should be considered
maxDelaySec = 1.25e-3;  

% Framing parameters
winSize = 2 * round(winSec * fs / 2);
hopSize = 2 * round(0.5 * winSec * fs / 2);
window  = hann(winSize);

% Relate maximum delay to samples (lags)                 
maxLag = ceil(maxDelaySec * fs);


%% 3. DECOMPOSE INPUT INTO INDIVIDUAL FREQUENCY CHANNELS
% 
% 
% Gammatone filtering    
bm_L = gammaFB(binaural(:,1),fs,fRange(1),fRange(2));
bm_R = gammaFB(binaural(:,2),fs,fRange(1),fRange(2));

% Number of auditory filters
nFilter = size(bm_L,2);

% Allocate memory
itd = zeros(nFilter,N);


%% 4. FRAME-BASED CROSS-CORRELATION ANALYSIS
% 
% 
% Loop over number of auditory filters
for ii = 1 : nFilter
    
    % Framing
    frames_L = frameData(bm_L(:,ii),winSize,hopSize,window,false);
    frames_R = frameData(bm_R(:,ii),winSize,hopSize,window,false);
    
    % Cross-correlation analysis to estimate ITD
    [CCF,lags] = xcorrNorm(frames_L,frames_R,maxLag);
    
    % Integrate cross-correlation pattern across all frames
    CCFsum = mean(CCF,2);

    % Estimate interaural time delay 
    itd(ii,:) = findITD(CCFsum,fs,lags,N);
end