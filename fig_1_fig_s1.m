% Estimate regression to the mean bias and regression dilution bias using
% synthetic observations. This script accompanies the paper
% Pitkanen et al., 2016, Artificial bias typically neglected in comparisons
% of uncertain atmospheric data, GRL
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

% create or load synthetic observations, this if clause is just a manual
% switch

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

% number of data points in x_true_original
x_true_n = 1e6;

% center and standard deviation of the gaussian distribution
x_true_center = 0.6;
x_true_std = 0.1;

% initialize random number generator
random_seed = 10;
rng(random_seed)

% values for x_true
x_true_original = randn(x_true_n,1) .* x_true_std + x_true_center;

% generate the synthetic observations.
% x_true is a cell variable containing x_true_original for each
% synthetic data set, which makes handling of variables
% consistent between x1, x2 and x_true
[ x_true, x1, dx1_abs, dx1_rel, x2, dx2_abs, dx2_rel ] = ...
    generate_observations( ...
        x_true_original, ...
        uncertainties1, dx1_bias, ...
        uncertainties2, dx2_bias);


%% plot data
% paper figure 1
if 1==1 
    
    % these bins are used for the solid lines representing ratio and
    % difference
    bin_edges   = [-inf, 0.15:.05:1.05, inf];
    bin_centers = mean([ bin_edges(1:(end-1)); bin_edges(2:end)]);
    
    % these are
    box_edges   = [-inf, 0.4, 0.6, 0.8, inf];
    box_centers = [ 0.3, 0.5, 0.7, 0.9];
 
    % eps becomes a bit smaller than grl full page recommendations
    set_fig_properties(12,20,8,1,4)
    fontsize_legend = get(0,'DefaultAxesFontSize')+2;
    fontsize_title  = get(0,'DefaultAxesFontSize');
        
    % this will position the subplot labels (a), (b), etc.
    subplot_label_pos = [0.04, 0.98];
        
    % choose which data sets to use
    i_set=3;

    % limits for the x axis
    x_lim       = [  0.0, 1.2];  
    
    % limits for the y-axis in ratio and difference plots
    y_lim_ratio = [ 0.3, 2.6];
    y_lim_diff  = [-0.4, 0.4];
  
    
    % histogram -----------------------------------------------------------
    n_rows = 3;
    n_cols = 2;
    figure
    subplot(n_rows, n_cols, 1)

    % x1
    [N1,~]=histc(x1{i_set}, bin_edges);
    h_x1 = stairs(bin_edges,N1,'r-','linewidth',4);
    hold on

    % x2
    [N2,~]=histc(x2{i_set}, bin_edges);
    h_x2 = stairs(bin_edges,N2,'b-','linewidth',2);
    
    % true data
    [Ntrue,~]=histc(x_true{i_set}, bin_edges);
    h_true = stairs(bin_edges,Ntrue, 'k-', 'linewidth',2);
    
    xlabel(x1_label{i_set})
    ylabel('Number of data points')
    
     
    legend([h_true, h_x1, h_x2], ...
           {'x_{true}', ...
           x1_label{i_set}, ...
           'x_{2}'}, ...
           'fontsize', fontsize_legend)


    grid on
    set(gca,'xlim',x_lim, 'ylim', [0, 3e5])
    axis square
    text(subplot_label_pos(1),subplot_label_pos(2),'(a)',...
         'Units', 'Normalized', ...
         'VerticalAlignment', 'Top', ...
         'fontsize',12)
    
    % scatter x1 vs x2 ----------------------------------------------------
    subplot(n_rows, n_cols, 2)
        
    [ber,mer] = bivariate_fit( x1{i_set}, x2{i_set}, ...
                               sqrt(0.05.^2 + (0.15 .* x1{i_set}).^2), ...
                               sqrt(0.05.^2 + (0.15 .* x2{i_set}).^2), 0, 0.79);
    a = polyfit(x1{i_set}, x2{i_set}, 1);

    % cloudplot
    cloudPlot(x1{i_set}, x2{i_set}, [x_lim, x_lim], 0, [100, 100])

    caxis([0,8e2])

    hold on
    plot(x_lim,x_lim,'k-','linewidth',0.5);
    
    
    % boxes
    [x1_binned, x2_binned] = bin_data(x1{i_set}, x2{i_set}, box_edges);
    
    boxPlot_cell(x2_binned, box_centers, [], 0.05, 2)
    hold on
    
    % fits
    ax_polyfit     = plot(bin_edges, a(1) .* bin_edges + a(2), 'b--');
    ax_bivar_fit   = plot(bin_edges, mer.*bin_edges + ber,'ko-');   
    
    % legend with fit parameters
    legend([ax_polyfit, ax_bivar_fit], ...
        {sprintf('OLS slope = %4.2f', a(1)) , ...
        sprintf('Bivariate slope = %4.2f',mer)}, ...
        'location','southeast', ...
        'fontsize', fontsize_legend);

    x1_str = sprintf('%s', x1_label{i_set}); 
    x2_str = sprintf('x_{2}'); 

    xlabel(x1_str)
    ylabel(x2_str)

    set(gca,'xlim',x_lim,'ylim',x_lim)
    grid on
    axis square
    text(subplot_label_pos(1),subplot_label_pos(2),'(b)',...
         'Units', 'Normalized', ...
         'VerticalAlignment', 'Top', ...
         'fontsize',30)
    
    % scatter x_true vs x2 ------------------------------------------------
    subplot(n_rows, n_cols, 3)
    
    axis square

    [h, canvas] = cloudPlot(x_true{i_set}, x2{i_set}, [x_lim, x_lim], 0, [100, 100]);
          
    caxis([0,2e3])
  
    hold on
    plot(x_lim, x_lim, 'k-', 'linewidth', 0.5);
    
    % plot line fits
    % calculate bivariate with proper uncertainties for x2, relative
    % uncertainty will change the slope a little
    [bo,mo] = bivariate_fit( x_true{i_set}, x2{i_set}, ...
                             1e-9.*ones(numel(x1{i_set}),1), ...
                             sqrt(0.05.^2 + (0.15 .* x2{i_set}).^2), 0, 1.0);
                         
    a = polyfit(x_true{i_set}, x2{i_set}, 1);

    ax_ols_fit   = plot(bin_edges, a(1) .* bin_edges + a(2), 'b--');
    hold on
    ax_bivar_fit = plot(bin_edges, mo.*bin_edges     + bo,   'ko-');

    legend([ax_ols_fit, ax_bivar_fit], ...
            {sprintf('OLS slope = %.2f', a(1)) , ...
            sprintf('Bivariate slope = %.2f',mo)}, ...
            'location','southeast', ...
            'fontsize', fontsize_legend);

    x1_str = 'x_{true}'; 
    x2_str = 'x_{2}'; 
    xlabel(x1_str)
    ylabel(x2_str)

    set(gca,'xlim',x_lim,'ylim',x_lim)
    grid on
    axis square
    text(subplot_label_pos(1),subplot_label_pos(2),'(c)',...
     'Units', 'Normalized', ...
     'VerticalAlignment', 'Top', ...
     'fontsize',12)
 
    
    % difference x2 / x1 --------------------------------------------------
    subplot(n_rows, n_cols, 4)
    
    % cloudplot
    cloudPlot(x1{i_set}, x2{i_set} - x1{i_set}, [x_lim, y_lim_diff], 0, [100, 100])
          
    caxis([0,5e2])

    hold on
    plot(x_lim, [0,0], 'k-', 'linewidth', 0.5)
    
    % boxes
    [~, diff_binned] = bin_data( x1{i_set}, x2{i_set} - x1{i_set}, box_edges);
    
    boxPlot_cell(diff_binned, box_centers,[],0.05,2)
    hold on
    
    % simulated lines
    [~, diff_binned] = bin_data(x1{i_set}, x2{i_set} - x1{i_set}, bin_edges);
    diff_binmedians = cellfun(@median, diff_binned);
    plot(bin_centers,diff_binmedians,'r')
    
    grid on
    ylabel('x_2 - x_1')
    xlabel(x1_label{i_set})
    
    set(gca,'xlim',x_lim,'ylim',y_lim_diff)
    axis square
    
    %text(subplot_label_pos(1),subplot_label_pos(2),'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
    text(subplot_label_pos(1),subplot_label_pos(2),'(d)',...
     'Units', 'Normalized', ...
     'VerticalAlignment', 'Top', ...
     'fontsize',12)
    
    % ratio x2 / x1 -------------------------------------------------------
    subplot(n_rows, n_cols, 5)
    
    cloudPlot(x1{i_set}, x2{i_set}   ./ x1{i_set}, ...
              [x_lim, y_lim_ratio], 0, [100, 100])
    

    caxis([0,1e3])

    hold on
    plot(x_lim, [1,1], 'k-', 'linewidth', 0.5)
    
    % boxes
    [~, ratio_binned] = bin_data(x1{i_set}, x2{i_set} ./ x1{i_set}, box_edges);
    
    boxPlot_cell(ratio_binned, box_centers,[],0.05,2)
    hold on
    
    % simulated line
    [~, ratio_binned] = bin_data(x1{i_set}, x2{i_set} ./ x1{i_set}, bin_edges);
    ratio_binmedians = cellfun(@median, ratio_binned);
    plot(bin_centers,ratio_binmedians,'r')
    
    grid on
    ylabel('x_{2} / x_{1}')
    xlabel(x1_label{i_set})
    
    set(gca,'xlim',x_lim,'ylim',y_lim_ratio)
    axis square
    text(subplot_label_pos(1),subplot_label_pos(2),'(e)',...
         'Units', 'Normalized', ...
         'VerticalAlignment', 'Top', ...
         'fontsize',12)
    
    % ratio x2 / x1--------------------------------------------------------
    subplot(n_rows, n_cols, 6)
    
    % calculate ratios x1-x2 for each data set
    ratio_binmedians = nan(numel(x1), numel(bin_centers));
    for l = 1:numel(x1)
        [~, ratio_binned]     = bin_data(x1{l}, x2{l} ./ x1{l}, bin_edges);
        ratio_binmedians(l,:) = cellfun(@median, ratio_binned);
    end
    
    plot(bin_centers, ratio_binmedians)
    hold on
    plot(x_lim, [1,1], 'k-', 'linewidth', 0.5)
    
    [hl, ho] = legend(x1_label, 'fontsize', fontsize_legend-4);
        
    xlabel('x_{1}')
    ylabel('x_{2} / x_{1}')
    grid on
    
    set(gca,'xlim',x_lim, 'ylim', y_lim_ratio)
    
    % add subplot label
    axis square
    
    text(subplot_label_pos(1),subplot_label_pos(2),'(f)',...
         'Units', 'Normalized', ...
         'VerticalAlignment', 'Top', ...
         'fontsize',12)     
    
    % change colormap, this applies for all subplots
    cmap = get_cloudPlot_canvas_percentile_cmap(canvas);

end


% paper figure S1
if 1==0 

    % set_fig_properties(width_cm,height_cm,fontsize_pt,linewidth,markersize)
    set_fig_properties(7,7,8,0.5,6)
    fontsize_legend = get(0,'DefaultAxesFontSize');
    
    i_set=3;
    
    % set appropriate limits for the data
    extreme_limits = [  0.2, 1.0];
    y_lim          = [ -0.2, 1.7];

    ind = x1{i_set} <= extreme_limits(1) | x1{i_set} >= extreme_limits(2);
    size_extreme = size( x_true{i_set}(ind));
    
    % values on x-axis where the data points go
    box_x = [1, 1.5, 2];
    
    % filter extreme subgroups
    ind1      = x1{i_set} <= extreme_limits(1) | x1{i_set} >= extreme_limits(2);
    x1_ex     = x1{i_set}(ind1);
    x2_ex     = x2{i_set}(ind1);
    x_true_ex = x_true{i_set}(ind1);
    
    % the function will draw a sample from the input and return the indices
    % for the sample
    ind_ex = randperm(numel(x_true_ex), 1000);

    % plot the connecting lines for the sampled points
    plot(box_x, [x1_ex(ind_ex), x_true_ex(ind_ex), x2_ex(ind_ex)], ' k-')
    hold on

    % plot the points
    plot(ones( size( x1{i_set}(ind))),          x1{i_set}(ind),    ' r.')
    plot(ones( size( x2{i_set}(ind))).*2.0,     x2{i_set}(ind),    ' b.')
    plot(ones( size( x_true{i_set}(ind))).*1.5, x_true{i_set}(ind),' k.')

    set(gca,'xlim',[1,2],'ylim',y_lim,'ygrid','on')
    
    x1_str = sprintf('%s', x1_label{i_set}); 
    x2_str = 'x2';
    x_true_str = 'x_{true}'; 
    set(gca,'xtick',[1,1.5,2],'xticklabel',{x1_str, x_true_str, x2_str})

    ylabel('Extreme subgroup value')

end
