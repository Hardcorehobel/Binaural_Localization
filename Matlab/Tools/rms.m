function y = rms(x, dim)
%RMS    Root mean squared value.
%   For vectors, RMS(X) is the root mean squared value in X. For matrices,
%   RMS(X) is a row vector containing the RMS value from each column. For
%   N-D arrays, RMS(X) operates along the first non-singleton dimension.
%
%   Y = RMS(X,DIM) operates along the dimension DIM.

if nargin==1
  y = sqrt(mean(x .* conj(x)));
else
  y = sqrt(mean(x .* conj(x), dim));
end
