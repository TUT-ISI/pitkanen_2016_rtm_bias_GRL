function [] = fig_s3()
% Creates an illustration of regression dilution for synthetic observations
% made from autocorrelating true values x_true. Also the effect of temporal
% averaging is shown. The figure is a reproduction of the Fig. S3 in the
% paper Pitkanen et al., 2016, Artificial bias typically neglected in
% comparisons of uncertain atmospheric data, GRL. This script has minimal
% error checking, so make sure you understand the code you run.
%
% INPUT:    none
%   
% OUTPUT:   none
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

clear all
close all

% the number of values in each sample. subsampling will be done
% n_sample_sets times for each n_samples value
n_samples =[1, 2, 3, 4, 5, 8, 11, 15, 20, 25, 30];
n_sample_sets = 10000;

% the number of values in the x_true
x_true_n = 1e6;

% autocorrelation of x_true with lag 1
roo1 = 0.0;
[ x1, x2 ] = generate_autocorrelated_observations( roo1, x_true_n );
[a_na, ~, ~] = subsample_and_fit(x1{1},x2{1},n_samples,n_sample_sets);


% autocorrelation of x_true with lag 1
roo2 = 0.5;
[ x1, x2 ] = generate_autocorrelated_observations( roo2, x_true_n );
[a_a50, ~, ~] = subsample_and_fit(x1{1},x2{1},n_samples,n_sample_sets);


% autocorrelation of x_true with lag 1
roo3 = 0.95;
[ x1, x2 ] = generate_autocorrelated_observations( roo3, x_true_n );
[a_a95, ~, ~] = subsample_and_fit(x1{1},x2{1},n_samples,n_sample_sets);


% -------------------------------------------------------------------------
% set the default figure properties
set_fig_properties(7.5,7.5,8,3,10)

figure

plot(n_samples, a_na(:,1), 'b.-')
hold on
plot(n_samples, a_a50(:,1), 'r.-')
plot(n_samples, a_a95(:,1), 'k.-')

ylabel('OLS slope')
xlabel('N data points per average')

leg_ax = legend({'Autoc. 0.00','Autoc. 0.50, lag 1','Autoc. 0.95, lag 1'}, ...
                'location','southeast');
leg_pos = get(leg_ax,'position');
set(leg_ax, 'position', [leg_pos(1).*0.9, leg_pos(2) .* 2, leg_pos(3), leg_pos(4)])
grid on
set(gca,'xlim', [0,30], 'ylim', [0.4, 1])

end


function [ x1, x2 ] = generate_autocorrelated_observations( roo, x_true_n )
%[ x1, x2 ] = generate_autocorrelated_observations( roo, x_true_n )
%
% Creates an autocorrelating x_true_original and then uses
% generate_observations() to create synthetic observations from
% x_true_original
%
% INPUT:
%   roo         autocorrelation of x_true_original, lag 1
%   x_true_n    number of values in x_true_original
%
% OUTPUT:
%   x1, x2   synthetic observations of x_true_original

% set the absolute and relative uncertainties for x1 corresponding
% some instruments or data products, etc.
% uncertainties1 = { ...
%   [sigma_abs, sigma_rel], ...
%   [sigma_abs, sigma_rel], ... etc}
% where sigma_abs is the absolute uncertainty and sigma_rel is the
% relative uncertainty
uncertainties1 = { [0.05, 0.15]};

% generate x2 with the same uncertainties (but a different
% instance of noise)
uncertainties2 = uncertainties1;

% set systematic biases for the data x1 and x2. Setting all to zero
% results in unbiased data
dx1_bias = {0};
dx2_bias = {0};

% -------------------------------------------------------------------------
% center and standard deviation of the gaussian distribution
x_true_center = 0.6;
x_true_std = 0.1;

% initialize random number generator
random_seed = 10;
rng(random_seed)

x_true_original = nan(x_true_n, 1);
x_true_original(1) = randn(1, 1) .* x_true_std;

% calculate an autocorrelating x_true_original
for i = 2:x_true_n
    x_true_original(i) = randn(1,1) + roo .* x_true_original(i-1);
end

% set standard deviation to x_true_std and set median value to
% x_true_center
x_true_original = x_true_original ./ std(x_true_original) .*...
                  x_true_std       + x_true_center;

% -------------------------------------------------------------------------
% generate the synthetic observations.
% x_true is a cell variable containing x_true_original for each
% synthetic data set, which makes handling of variables
% consistent between x1, x2 and x_true
[ ~, x1, ~, ~, x2, ~, ~ ] = ...
    generate_observations( ...
        x_true_original, ...
        uncertainties1, dx1_bias, ...
        uncertainties2, dx2_bias);
end




function [a,b,m] = subsample_and_fit(x1,x2,n_samples,n_sample_sets)
% calculates ordinary least squares (ols) and bivariate linear fit
% parameters for temporally averaged synthetic observations x1 and x2
%
% INPUT
%   x1              a synthetic observation created in main_rtm.m
%   x2              another synthetic observation created in main_rtm.m
%   n_samples       number of x1 points in each sample
%   n_sample_sets   number of samples
%
% OUTPUT
%   a               [OLS_slope, OLS_y_intercept]
%   b               bivariate slope
%   m               bivariate y-intercept

% initialize
a = nan(numel(n_samples),2);
b = nan(numel(n_samples),1);
m = nan(numel(n_samples),1);

% loop the number of data points in each subsample
for k = 1:numel(n_samples)

    % initialize
    samples_x1     = nan(n_sample_sets, n_samples(k));
    samples_x2     = nan(n_sample_sets, n_samples(k));

    % loop samples
    for l = 1:n_sample_sets

        % sample indices for consecutive values
        ind = (1:n_samples(k)) + n_samples(k) .* (l-1);

        % pick the sample
        samples_x1(l,:)     = x1(ind);
        samples_x2(l,:)     = x2(ind);
    end

    % calculate subsample means
    sample_x1_medians     = mean(samples_x1,2);
    sample_x2_medians     = mean(samples_x2,2);

    % make linear fits to subsample means. these are not needed for the
    % figure, but you can check that your bivariate fit always
    % gives a slope close to 1. The uncertainty estimates here are the same
    % as for unaveraged x1 and x2. This may not be entirely true, but for
    % this purpose its mostly important that they are the same for the
    % averaged x1 and x2.
    [b(k), m(k)] = bivariate_fit( ...
        sample_x1_medians, sample_x2_medians, ...
        sqrt(0.05.^2 + (0.15 .* sample_x1_medians).^2), ...
        sqrt(0.05.^2 + (0.15 .* sample_x2_medians).^2), 0, 1);

    a(k,:) = polyfit(sample_x1_medians, sample_x2_medians, 1);
end

end

