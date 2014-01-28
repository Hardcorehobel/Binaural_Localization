function [bin,hrtf] = filterHRTF(in,fs,azim,select)
%filterHRTF   Spatialize input signal using HRTF filtering.
%
%USAGE
%             BIN = filterHRTF(IN,FS);
%      [BIN,HRTF] = filterHRTF(IN,FS,AZIM,SELECT);
%
%INPUT PARAMETERS
%       IN : input signal [nSamples x 1]
%       FS : sampling frequency in Hertz
%     AZIM : azimuth angle [-90:5:90]
%   SELECT : string specifying HRTF catalog
% 
%            SURREY database using the Head And Torso Simulator (HATS)
%            [Hummersone, IEEE TASLP 2010]
% 
%            'SURREY_A'      = Anechoic 1.4m   197 nPoints
%            'SURREY_ROOM_A' = T60 = 0.32s    6259 nPoints
%            'SURREY_ROOM_B' = T60 = 0.47s   16001 nPoints
%            'SURREY_ROOM_C' = T60 = 0.68s   16001 nPoints
%            'SURREY_ROOM_D' = T60 = 0.89s   16349 nPoints
%
%OUTPUT PARAMETERS
%      BIN : binaural audio signal 
%     HRTF : head-related transfer functions [nPoints x 2]
%
%NOTE 
%   For proper functioning, the HRTF database is assumed to be located
%   in the subfolder HRTF_SURREY, relative to the position of this
%   function.  
%
%EXAMPLE
%   % Create white noise 
%   in = randn(32e3,1);
% 
%   % Filter signal with HRTFs corresponding to 45 degrees azimuth
%   out = filterHRTF(in,16E3,45);

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/01/21
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
if nargin < 4 || isempty(select); select = 'SURREY_A'; end


%% 2. LOAD HRTFs  
% 
% 
% Read HRTFs 
switch upper(select)
    case 'SURREY_A'
        % Root directory
        rootHRTF = [pwd,filesep,'Tools',filesep,'HRTF_SURREY',filesep,...
                    'Anechoic',filesep,'16kHz',filesep];
        % HRTF Name
        nameHRTF = ['CortexBRIR_0s_',num2str(azim),'deg_16k.wav'];
        
        % HATS database 
        [hrtf,fsRef] = wavread([rootHRTF,nameHRTF]);
        
    case 'SURREY_ROOM_A'
        % Root directory
        rootHRTF = [pwd,filesep,'Tools',filesep,'HRTF_SURREY',filesep,...
                    'Room_A',filesep,'16kHz',filesep];
        % HRTF Name
        nameHRTF = ['CortexBRIR_0_32s_',num2str(azim),'deg_16k.wav'];
        
        % HATS database
        [hrtf,fsRef] = wavread([rootHRTF,nameHRTF]);
        
    case 'SURREY_ROOM_B'
        % Root directory
        rootHRTF = [pwd,filesep,'Tools',filesep,'HRTF_SURREY',filesep,...
                    'Room_B',filesep,'16kHz',filesep];
        % HRTF Name
        nameHRTF = ['CortexBRIR_0_47s_',num2str(azim),'deg_16k.wav'];
        
        % HATS database
        [hrtf,fsRef] = wavread([rootHRTF,nameHRTF]);
        
    case 'SURREY_ROOM_C'
        % Root directory
        rootHRTF = [pwd,filesep,'Tools',filesep,'HRTF_SURREY',filesep,...
                    'Room_C',filesep,'16kHz',filesep];
        % HRTF Name
        nameHRTF = ['CortexBRIR_0_68s_',num2str(azim),'deg_16k.wav'];
        
        % HATS database
        [hrtf,fsRef] = wavread([rootHRTF,nameHRTF]);
        
    case 'SURREY_ROOM_D'
        % Root directory
        rootHRTF = [pwd,filesep,'Tools',filesep,'HRTF_SURREY',filesep,...
                    'Room_D',filesep,'16kHz',filesep];
        % HRTF Name
        nameHRTF = ['CortexBRIR_0_89s_',num2str(azim),'deg_16k.wav'];
        
        % HATS database
        [hrtf,fsRef] = wavread([rootHRTF,nameHRTF]);    
        
    otherwise
        error(['HRTF selection ''',upper(select),''' is not supported.'])
end

% HRTF Resampling
if ~isequal(fs,fsRef)
    hrtf = resample(hrtf, fs, fsRef);
end


%% 3. APPLY HRTFs
% 
% 
% Number of HRTF channels
nChanHRTF = size(hrtf,2);

% Allocate memory for binaural signal
bin = zeros(size(in,1),nChanHRTF);

% Filter input signal with HRTFs
for ii = 1 : nChanHRTF
    bin(:,ii) = fftfilt(hrtf(:,ii),in);
end


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************