function [itd,lags] = estimate_ITD_Broadband(binaural,fs,winSec,N)
%estimate_ITD_Broadband   Broadband ITD estimation.
%
% The interaural time difference (ITD) is estimated from binaural signals
% by computing the cross-correlation function (CCF) over short time frames
% (winSec). The resulting 2-dimensional cross-correlation function, which
% is a function of the number of time lags (ITDs) and the number of frames
% (nLags x nFrames), is integrated across all time frames and the most
% prominent peaks are assumed to reflect the estimated ITDs for N 
% sources.         
% 
%USAGE
%    [itd,lags] = estimate_ITD_Broadband(binaural,fs)
%    [itd,lags] = estimate_ITD_Broadband(binaural,fs,winSec,N)
%
%INPUT PARAMETERS
%   binaural : binaural input signal [nSamples x 2]
%         fs : sampling frequency in Hertz
%     winSec : frame size in seconds of the cross-correlation analysis
%              (default, winSec = 20E-3)
%          N : number of sources that should be localized (default, N = 1) 
%
%OUTPUT PARAMETERS
%        itd : estimated ITDs of all N sources [N x 1]
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
if nargin < 2 || nargin > 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameters
if nargin < 3 || isempty(winSec);   winSec   = 20E-3; end
if nargin < 4 || isempty(N);        N        = 1;     end


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


%% 3. SEGMENT INPUT SIGNAL INTO OVERLAPPING TIME FRAMES
% 
%  
% Framing
frames_L = frameData(binaural(:,1),winSize,hopSize,window,false);
frames_R = frameData(binaural(:,2),winSize,hopSize,window,false);


%% 4. CROSS-CORRELATION ANALYSIS
% 
% 
% Cross-correlation analysis 
[CCF,lags] = xcorrNorm(frames_L,frames_R,maxLag);

% Integrate cross-correlation pattern across all frames
CCFsum = mean(CCF,2);


%% 5. FIND ITD
% 
% 
% Find ITDs
itd = findITD(CCFsum,fs,lags,N);