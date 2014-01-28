function itdEst = findITD(Csum,fs,lags,nSources)


%% 1. CHECK INPUT ARGUMENTS  
% 
% 
% Check for proper input arguments
if nargin < 3 || nargin > 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 4 || isempty(nSources); nSources = 1; end


%% 2. FIND LOCAL PEAKS
% 
% 
% Find peaks, also consider endpoints as peak candidates
pIdx = findpeaks([0; Csum(:); 0]);
pIdx = pIdx - 1;

% Rank peaks
[temp,idx] = sort(Csum(pIdx),'descend'); %#ok

% Restrict number of ITD estimates
nEst = min(numel(idx),nSources);

% Integer lag: Take most salient peaks
lagInt = lags(pIdx(idx(1:nEst)));


%% 3. PERFORM INTERPOLATION
% 
% 
% Fractional lag: Refine peak position by parabolic interpolation
lagDelta = interpolateParabolic(Csum,pIdx(idx(1:nEst)));

% Stepsize of lags (required for interpolation)
deltaT = abs(diff(lags(1:2)));

% Final interaural time delay estimates
itdEst = (lagInt + (deltaT * lagDelta))/fs;    