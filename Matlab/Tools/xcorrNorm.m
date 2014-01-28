function [xcorr,lags] = xcorrNorm(inL,inR,maxLag,bDetrend,bNorm)
%xcorrNorm   Normalized time-domain cross-correlation function.
%
%USAGE
%   [XCORR,LAGS] = xcorrNorm(INL,INR)
%   [XCORR,LAGS] = xcorrNorm(INL,INR,MAXLAG,bDETREND,bNORM)
%   
%INPUT ARGUMENTS
%        INL : left input arranged as  [nSamples x nChannels]
%        INR : right input arranged as [nSamples x nChannels]
%     MAXLAG : computation is performned over the lag range -MAXLAG:MAXLAG
%             (default, MAXLAG = nSamples-1)
%   bDETREND : substract mean     (default, bDETREND = true)
%      bNORM : normalization flag (default, bNORM    = true)
%
%OUTPUT ARGUMENTS
%     XCORR : cross-correlation function [nSamples x nChannels] 
%      LAGS : time lags of cross-correlation function [2*MAXLAG+1 x 1]
% 
%REFERENCES
%   [1]  Roman, N., Wang, D. L. and Brown, G. J., "Speech segregation based
%        on sound localization", J. Acoust. Soc. Amer., vol. 114, no. 4,
%        pp. 2236-2252, 2003.
% 
%   See also acorrNorm and calcAMDF.   

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2009 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
% 
%   History :  
%   v.0.1   2009/10/10
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(maxLag);   maxLag   = size(inL,1)-1; end
if nargin < 4 || isempty(bDetrend); bDetrend = true;          end
if nargin < 5 || isempty(bNorm);    bNorm    = true;          end
    

%% **************************  MEX PROCESSING  ****************************
% 
% 
[xcorr,lags] = xcorrNormMEX(inL,inR,maxLag,bDetrend,bNorm);

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