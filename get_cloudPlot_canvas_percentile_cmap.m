function [ cmap ] = get_cloudPlot_canvas_percentile_cmap( canvas )
%get_cloudPlot_canvas_percentile_cmap 
%   When you make a cloud plot, this function will attempt to make the
%   shades of the colormap match with percentiles in the data density. If
%   you use this on a figure with several cloudplot in different subplots,
%   the same colomap will be used for all of them, and the percentiles
%   might not match with the shades.
% 
% Copyright (C) 2016 Mikko Pitkanen
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
%
% Input: canvas         an output parameter from cloudPlot.m

% number of data on the canvas
N = sum(canvas(:));

percentiles = [0.10, 0.25, 0.50, 0.75, 0.90];
N_percentiles = percentiles .* N;

% minimum density of data points on canvas
n_min = min(canvas(:));
% maximum density of data points on canvas
n_max = max(canvas(:));

n_steps = (n_max - n_min) / 100;

n_canvas_percentiles = nan(1,numel(percentiles));

% loop percentiles
for i = 1:numel(percentiles)
    
    n_current = n_max;
    N_current = 0;
    
    % find the right n for percentiles(i)
    while N_current < N_percentiles(i) 
        
        n_current = n_current - n_steps;
        
        ind = canvas >= n_current;
        
        N_current = sum(canvas(ind));
    
    end
    n_canvas_percentiles(i) = n_current;
    
end

caxis([0,n_max])

% four shades of gray
greys = [0.2;
         0.4;
         0.6;
         0.8];
          
% create the colormap
cmap = nan(1000,3);

cmap_percentiles = n_canvas_percentiles / n_max * 1000;

% below percentiles(1)
cmap(                         1:round(cmap_percentiles(5)), :) = 1;
% between percentiles(1) and percentiles(2)
cmap(round(cmap_percentiles(5)):round(cmap_percentiles(4)), :) = greys(4);
% between percentiles(2) and percentiles(3)
cmap(round(cmap_percentiles(4)):round(cmap_percentiles(3)), :) = greys(3);
% between percentiles(3) and percentiles(4)
cmap(round(cmap_percentiles(3)):round(cmap_percentiles(2)), :) = greys(2);
% between percentiles(4) and percentiles(end)
cmap(round(cmap_percentiles(2)):round(cmap_percentiles(1)), :) = greys(1);
% above percentiles(end)
cmap(round(cmap_percentiles(1)):end, :) = 0;
              
colormap(cmap)
end

