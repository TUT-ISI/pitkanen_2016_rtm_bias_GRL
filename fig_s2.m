% Estimate regression to the mean bias and regression dilution bias using
% synthetic observations. This script accompanies the paper
% Pitkanen et al., 2016, Artificial bias typically neglected in comparisons
% of uncertain atmospheric data, GRL. There is only minimal error checking
% so make sure you understand the code you run.
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

%% load or create the data in the article 

clear all
close all

% the number of labels, uncertainties and biases must be the same
% for x1 and x2
x1_label = { ...
    'x_{1}: Demo 1', ...
    'x_{1}: AERONET', ...
    'x_{1}: MODIS', ...
    'x_{1}: MISR', ...
    'x_{1}: SL501 lo cos.', ...
    'x_{1}: SL501 hi cos.', ...
    'x_{1}: Demo 2'};

% for simplicity x2 data sets are not given labels as specific as
% for x1
x2_label = { ...
    'x_{2}', ...
    'x_{2}', ...
    'x_{2}', ...
    'x_{2}', ...
    'x_{2}', ...
    'x_{2}', ...
    'x_{2}'};

% set the absolute and relative uncertainties for x1 corresponding
% some instruments or data products, etc.
% uncertainties1 = { ...
%   [sigma_abs, sigma_rel], ...
%   [sigma_abs, sigma_rel], ... etc}
% where sigma_abs is the absolute uncertainty and sigma_rel is the
% relative uncertainty
uncertainties1 = { ...
    [0.3, 0.0], ...     Demonstrational1
    [0.01, 0.0], ...    AERONET
    [0.05, 0.15], ...   MODIS
    [0.05, 0.20], ...   MISR
    [0.0, 0.07], ...    SL501
    [0.0, 0.16], ...    SL501
    [0.0, 0.3]}; %      Demonstrational2


% generate x2 with the same uncertainties (but a different
% instance of noise)
uncertainties2 = uncertainties1;

% set systematic biases for the data x1 and x2. Setting all to zero
% results in unbiased data
dx1_bias = {0,0,0,0,0,0,0};
dx2_bias = {0,0,0,0,0,0,0};


% load the predefined variable x_true_original. 
% it contains 1e6 values of randomly sampled AOT_500
% from lev 2.0 direct sun product for the stations
% GSFC, Lille, BONDVILLE and Kuopio. The data was downloaded
% from http://aeronet.gsfc.nasa.gov/ on April 26, 2016.

% We thank the PI investigators and their staff for
% establishing and maintaining the 4 sites used in this
% investigation.
load('./x_true_original_aod');

% number of data points in x_true_original
x_true_n = numel(x_true_original);

% generate the synthetic observations.
% x_true is a cell variable containing x_true_original for each
% synthetic data set, which makes handling of variables
% consistent between x1, x2 and x_true
[ x_true, x1, dx1_abs, dx1_rel, x2, dx2_abs, dx2_rel ] = ...
    generate_observations_fig_s2( ...
        x_true_original, ...
        uncertainties1, dx1_bias, ...
        uncertainties2, dx2_bias);


%% plot data
%  paper figure S2

% MODIS = 3
i = 3;

% MISR = 4
k = 4;

% these bins are used for the solid lines representing ratio and
% difference
bin_edges   = [-inf, 0:.05:1, inf];
bin_centers = mean([ bin_edges(1:(end-1)); bin_edges(2:end)]);

% these are used for binning the data with regard to the boxes
box_edges   = [-inf, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, inf];
box_centers = [-0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3];

% this will position the subplot labels (a), (b), etc.
subplot_label_pos = [0.04, 0.98];

cmap_log = false;

% eps becomes a bit smaller than grl full page recommendations
set_fig_properties(7,15,8,1.0,3)
fontsize_legend = get(0,'DefaultAxesFontSize');
fontsize_title  = get(0,'DefaultAxesFontSize');

x_lim      = [ 0.0, 1.2];

y_lim_diff = [-0.8, 0.8];

% ---------------------------------------------------------------------
subplot(2,1,1)

[ber,mer] = ...
    bivariate_fit( x1{i}, x1{k}, ...
                   sqrt(0.05.^2 + (0.15 .* x1{2}).^2), ...
                   sqrt(0.05.^2 + (0.20 .* x1{2}).^2), 0, 1.0);

a = polyfit(x1{i}', x1{k}', 1);

% scatter
[h, canvas] = cloudPlot(x1{i}, x1{k}, ...
                       [x_lim, x_lim], cmap_log, [100, 100]);

hold on
plot(x_lim,x_lim,'k-','linewidth',1);


% plot boxes
[~, x2_binned] = bin_data(x1{i}, x1{k}, box_edges);

boxPlot_cell(x2_binned, box_centers, [], 0.05, 1)
hold on

% plot fits
ax_polyfit     = plot(bin_edges, a(1) .* bin_edges + a(2), 'b--');
ax_bivar_fit   = plot(bin_edges, mer.*bin_edges + ber,'ko-');

legend([ax_polyfit, ax_bivar_fit], ...
        {sprintf('OLS slope = %4.2f', a(1)) , ...
         sprintf('Bivariate slope = %4.2f',mer)}, ...
        'location','southeast');

x1_str = sprintf('%s', x1_label{i});
x2_str = sprintf('%s', x1_label{k});

xlabel(x1_str)
ylabel(x2_str)

set(gca,'xlim',x_lim,'ylim',x_lim)
grid on

text(subplot_label_pos(1),subplot_label_pos(2),'(a)', ...
    'Units', 'Normalized', 'VerticalAlignment', 'Top')


% ---------------------------------------------------------------------
subplot(2,1,2)

% scatter cloudplot      
cloudPlot(x1{i}, x1{k} - x1{i}, [x_lim, y_lim_diff], cmap_log, [100, 100])

caxis([0,1e3])

hold on
plot(x_lim, [0,0], 'k-', 'linewidth', 0.5)

% boxes
[~, diff_binned] = bin_data(x1{i}, x1{k} - x1{i}, box_edges);

boxPlot_cell(diff_binned, box_centers,[],0.05,1)
hold on

% simulated lines
h_sim = plot(box_centers, cellfun(@median,diff_binned),'r');

[hl, ho] = legend(h_sim,{'Median difference'},'location','north');

grid on
ylabel([x1_label{k} ' - ' x1_label{i}])
xlabel(x1_label{i})

set(gca,'xlim',x_lim,'ylim',y_lim_diff)

text(subplot_label_pos(1),subplot_label_pos(2),'(b)', ...
    'Units', 'Normalized', 'VerticalAlignment', 'Top')

% change colormap, this applies for all subplots
cmap = get_cloudPlot_canvas_percentile_cmap(canvas);  

