function azimEst = findAzimuth(CWarped,azimuth,nSources)


%% 1. CHECK INPUT ARGUMENTS  
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 3 || isempty(nSources); nSources = 1; end


%% 2. FIND LOCAL PEAKS
% 
% 
% Find peaks, also consider endpoints as peak candidates
pIdx = findpeaks([0; CWarped(:); 0]);
pIdx = pIdx - 1;

% Rank peaks
[temp,idx] = sort(CWarped(pIdx),'descend'); %#ok

% Number of azimuth estimates
nEst = min(numel(idx),nSources);

% Integer azimuth: Take most significant peaks
azimInt = azimuth(pIdx(idx(1:nEst)));


%% 3. PERFORM INTERPOLATION
% 
% 
% Fractional azimuth: Refine peak position by parabolic interpolation
azimDelta = interpolateParabolic(CWarped,pIdx(idx(1:nEst)));

% Azimuth stepsize (required for interpolation)
deltaT = abs(diff(azimuth(1:2)));

% Final interaural time delay estimates
azimEst = azimInt(:) + (deltaT * azimDelta(:));