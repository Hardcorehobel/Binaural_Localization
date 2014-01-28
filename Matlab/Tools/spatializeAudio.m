function binaural = spatializeAudio(audio,fs,azimuth,room)

% Check for proper input arguments
if nargin ~= 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Number of audio files
[nSamples,nSources] = size(audio);

% Allocate memory
binaural = zeros(nSamples,2);

% Loop over number of audio files
for ii = 1 : nSources
    % Spatialize signal using HRTF processing
    binaural = binaural + filterHRTF(audio(:,ii),fs,azimuth(ii),room);
end