function [ x_true, x1, dx1_abs, dx1_rel, ...
                   x2, dx2_abs, dx2_rel ] = ...
     generate_observations( x_true_original, ...
                            uncertainties1, dx1_bias, ...
                            uncertainties2, dx2_bias)

%FUNCTION generate_observations(x_true_original, ...
%                               uncertainties1, dx1_bias, ...
%                               uncertainties2, dx2_bias)
%
% Returns synthetic observations made from a random number data set. See
% Pitkanen et al., 2016, GRL for additional explanations. This script has
% minimal error checking, so make sure you understand the code you run.
%
%INPUT:
%   x_true_original     values of x_true, has to be a vector 
%   uncertainties1      for x1: {[sigma_abs, sigma_rel],[sigma_abs, sigma_rel], ...}
%                       each data set x1 requires an uncertainty [sigma_abs, sigma_rel]
%   uncertainties2      same as uncertainties1 but for x2
%   dx1_bias            systematic bias for synthetic observation x1
%   dx2_bias            systematic bias for synthetic observation x2
%
%OUTPUT
%   x_true              cell variable that contains each set of x_true
%   x1                  synthetic observations x1 of each set of x_true
%   dx1_abs             difference x1 - x_true caused by uncertainty1 sigma_abs
%   dx1_rel             difference x1 - x_true caused by uncertainty1 sigma_rel
%   x2                  synthetic observations x2 of each set of x_true
%   dx2_abs             difference x2 - x_true caused by uncertainty2 sigma_abs
%   dx2_rel             difference x2 - x_true caused by uncertainty2 sigma_rel
%
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

%% check inputs
if     ~( isvector(x_true_original) )
    error('Check input x_true_center in function generate_observations()')

elseif ~( iscell(uncertainties1) && ...
          iscell(dx1_bias)       && ...
          iscell(uncertainties2) && ...
          iscell(dx2_bias))
    error('Check input type in function generate_observations()')
    
elseif ~( numel(uncertainties1) == numel(dx1_bias)       && ...
          numel(uncertainties1) == numel(uncertainties2) && ...
          numel(uncertainties1) == numel(dx2_bias))
    error('Check input dimensions in function generate_observations()')
      
end


%% generate the data

% initialize variables
number_of_data_sets = numel(uncertainties1);
x_true              = cell( 1, number_of_data_sets);
x1                  = cell( 1, number_of_data_sets);
dx1_abs             = cell( 1, number_of_data_sets);
dx1_rel             = cell( 1, number_of_data_sets);
x2                  = cell( 1, number_of_data_sets);
dx2_abs             = cell( 1, number_of_data_sets);
dx2_rel             = cell( 1, number_of_data_sets);

% loop through each set of x1
for i = 1:number_of_data_sets
    
    % use the same x_true for each set of x1 and x2
    x_true{i}   = x_true_original;
    
    % generate the synthetic data x1
    [x1{i}, dx1_abs{i}, dx1_rel{i}] = ...
        create_synthetic_data_set(x_true{i}, uncertainties1{i}(1), uncertainties1{i}(2));

    % generate the synthetic data x2
    [x2{i}, dx2_abs{i}, dx2_rel{i}] = ...
        create_synthetic_data_set(x_true{i}, uncertainties2{i}(1), uncertainties2{i}(2));
    
    % add systematic bias to data
    x1{i} = x1{i} + dx1_bias{i};
    x2{i} = x2{i} + dx2_bias{i};

end

end


function [ x, dx_abs, dx_rel] = create_synthetic_data_set(x_true, uncert_abs, uncert_rel)
% [ x, dx_abs, dx_rel] = create_synthetic_data_set(x_true, uncert_abs, uncert_rel)
%    Creates a synthetic observation data set based on x_true by adding
%    absolute and relative observation uncertainty uncert_abs and
%    uncert_rel. uncert_rel can be defined with regard to x_true or x, so
%    remember to choose either one!
%
% Input:
%   x_true      true values representing an error free physical quantity
%   uncert_abs  absolute random uncertainty
%   uncert_rel  relative random uncertainty
%
% Output:
%   x           a synthetic observation of x
%   dx_abs      gaussian absolute random error that was added
%   dx_rel      gaussian relative random error that was added
%

% dimensions of x_true
dim_true = size(x_true);

% rand_dx_rel is the (0,1) normally distributed random number set for the
% relative uncertainty. rand_dx_abs is the same for the absolute
% uncertainty
rand_dx_abs = randn(dim_true);
rand_dx_rel = randn(dim_true);

% dx_abs is easy
dx_abs = rand_dx_abs .* uncert_abs;

% this if clause is a manual switch for choosing between different
% definitions of the relative uncertainty sigma_rel

% define relative error with regard to x
if 1==0
    % 1:   x       = x_true + dx_abs + dx_rel + dx_bias
    % 2:   dx_abs  = rand_dx_abs * uncert_abs
    % 3:   dx_rel  = rand_dx_rel * uncert_rel * x
    % 4:   dx_bias = 0
    % Solve x:
    x = ( dx_abs + x_true) ./ (1 - rand_dx_rel .* uncert_rel );

    dx_rel = rand_dx_rel .* uncert_rel .* x;
    
    
% define relative error with regard to x_true. dx_abs and dx_rel can be
% positive and negative and so they can cancel each others sometimes, which
% will reduce the variance of the synthetic observation. this was used in
% the paper Pitkanen et al., 2016, GRL
elseif 1==1
    
    % 1:   x       = x_true + dx_abs + dx_rel + dx_bias
    % 2:   dx_abs  = rand_dx_abs * uncert_abs
    % 3:   dx_rel  = rand_dx_rel * uncert_rel * x_true
    % 4:   dx_bias = 0

    dx_rel = rand_dx_rel .* uncert_rel .* x_true;
    x      = x_true + dx_abs + dx_rel;
        
% this definition of total uncertainty will create MODIS-like (Levy et
% al, 2013) synthetic observations. the absolute and random uncertainty do
% not cancel each others, which results in larger total uncertainty than in
% the previous option. the total uncertainty here is (uncert_abs +
% uncert_rel.*x_true) while in the previous the combined total uncertainty
% is sqrt(uncert_abs^2 + (uncert_rel.*x_true)^2). Compare your synthetic
% AERONET AOD and synthetic MODIS AOD with the results in Levy et al, 2013,
% fig 11 and you should get similar amount of data in the error envelope.
% This option was also used in the paper Pitkanen et al., 2016, GRL
elseif 1==1
    
    % 1:   x         = x_true + dx_random + dx_bias
    % 2:   dx_random = rand_dx * (uncert_abs + uncert_rel.*x_true)
    % 3:   dx_bias   = 0

    rand_dx   = randn(dim_true);
    dx_random = rand_dx .* (uncert_abs + uncert_rel .* x_true);
    dx_bias   = 0;
    x         = x_true + dx_random + dx_bias;
    
    dx_abs = nan(dim_true);
    dx_rel = nan(dim_true);
end

end



