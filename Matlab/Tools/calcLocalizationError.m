function [locErrDeg,newIdx] = calcLocalizationError(azRef,azEst)
%calcLocalizationError   Compute absolute localization error.
% 
%   For each true sound source position, the best matching localization
%   estimate is found which minimizes the absolute error.
%
%USAGE 
%   [locErrDeg,newIdx] = calcLocalizationError(azRef,azEst)
%
%INPUT ARGUMENTS
%   azRef : vector with true sound source positions      [K x 1]
%   azEst : vector with estimated sound source positions [N x 1]
%
%OUTPUT ARGUMENTS
%   locErrDeg : absolute localization error [N x 1]
%      newIdx : new indices 

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/01/26
%   ***********************************************************************

% Check for proper input arguments
if nargin ~= 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Number of estimated and reference sound source positions
nSourcesEst = numel(azEst);
nSourcesRef = numel(azRef);

% Allocate memory
estIdx    = 1:nSourcesEst;
refIdx    = 1:nSourcesRef;
newIdx    = zeros(nSourcesEst,1);
locErrDeg = zeros(nSourcesEst,1);

% Loop over all localization estimates and iteratively select the one with
% the minimum deviation from the reference position. 
for ii = 1:nSourcesEst
    
    currEst = azEst(estIdx);
    currRef = azRef(refIdx);
        
    currEstP = repmat(currEst(:).',[length(currRef) 1]);
    currRefP = repmat(currRef(:),  [1 length(currEst)]);
    
    estPIdx = repmat((1:length(estIdx)),  [length(currRef) 1]);
    refPIdx = repmat((1:length(refIdx)).',[length(currEst) 1]);
        
    % Compute localization error for all possible combinations
    currCosts = abs(currRefP(:) - currEstP(:));
    
    % Find best matching localization estimate
    [minErr,minIdx] = min(currCosts(:));
    
    % Store localization error
    locErrDeg(estIdx(estPIdx(minIdx))) = minErr;
    
    % Store new index 
    newIdx(ii) = estIdx(estPIdx(minIdx));
    
    % Remove index 
    estIdx(estPIdx(minIdx)) = [];
    refIdx(refPIdx(minIdx)) = [];
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