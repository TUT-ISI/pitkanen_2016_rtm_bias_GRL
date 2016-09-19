function [ x_binned, y_binned ] = bin_data( x, y, bin_edges )
% bin_data( x, y, bin_edges )
%
% This function bins x and y according to x. 
% 
% Input:
%   x           vector or matrix. reference data to use for binning
%   y           vector or matrix, same size as x. y will be binned
%               according to corresponding data points in x 
%   bin_edges   edges for bins
%
% Output:
%   x_binned    cell variable containing x data in bins bin_edges
%   y_binned    cell variable containing y data for which x is in bins
%               bin_edges
%
% % Copyright (C) 2016 Mikko Pitkanen
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check input
if ~(isnumeric(x) && isnumeric(y) && isnumeric(bin_edges))
    error('All input must be numeric in bin_data( x, y, bin_edges )')
    
elseif ~( isequal( size(x), size(y)))
    error('x and y must have the same size in bin_data( x, y, bin_edges )')
    
elseif ~isvector( bin_edges) && (numel(bin_edges) < 2)
    error('bin_edges must a vector with at least 2 elements in bin_data( x, y, bin_edges )')
    
end
    
% initialize
x_binned = cell(1, numel(bin_edges) - 1);
y_binned = cell(1, numel(bin_edges) - 1);

% loop bins
for i = 1:(numel(bin_edges) - 1)
    
    % find indices for x data in the current bin
    ind = x >= bin_edges(i) & ...
          x <  bin_edges(i+1);
    
    % place data into the bins, if any data was found
    if any(ind)
        x_binned{i} = x(ind);
        y_binned{i} = y(ind);
    else
        x_binned{i} = nan;
        y_binned{i} = nan;
    end
    
end

end

