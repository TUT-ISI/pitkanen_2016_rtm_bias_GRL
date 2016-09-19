function [a, b, S] = bivariate_fit( xi, yi, dxi, dyi, ri, b0 )
% Makes a linear bivariate fit to xi, yi data using York et al. (2004)
% algorithm:
%
% York, D et al., Unified equations for the slope, intercept, and standard
% errors of the best straight line, American Journal of Physics, 2004, 72,
% 3, 367-375, doi = 10.1119/1.1632486
%
% See especially Section III and Table I. The enumerated steps below are
% citations to Section III
% 
% Input variables:
%   xi, yi      x and y data points
%   wxi, wyi    errors for the data points xi, yi
%   ri          correlation coefficient for the weights
%   m0          initial guess m
%
% Output variables:
%   a           y-intercept, y = a +bx
%   b           slope
%
% Usage:
% [a, b] = bivariate_fit( xi, yi, dxi, dyi, ri, b0 )
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


% (1) Choose an approximate initial value of b
b = b0;

% (2) Determine the weights wxi, wyi, for each point.
wxi = 1 ./ dxi.^2;
wyi = 1 ./ dyi.^2;

alphai = sqrt(wxi .* wyi);
b_diff = 999;

% tolerance for the fit, when m changes by less than tol for two
% consecutive iterations, fit is considered found
tol = 1e-8;

% iterate until m changes less than tol
while abs(b_diff) >= tol

    b_prev = b;

    % (3) Use these weights wxi, wyi to evaluate Wi for each point.  
    Wi  = (wxi .* wyi) ./ (wxi    + b.^2 .* wyi - 2*b*ri*alphai);

    % (4) Use the observed points (xi ,yi) and Wi to calculate x_bar and
    % y_bar, from which Ui and Vi , and hence betai can be evaluated for
    % each point
    x_bar = sum(Wi .* xi) ./ sum(Wi);
    y_bar = sum(Wi .* yi) ./ sum(Wi);

    Ui = xi - x_bar;
    Vi = yi - y_bar;

    betai = Wi .* (Ui ./ wyi + b.* Vi ./ wxi - (b.*Ui + Vi) .* ri ./ alphai);

    % (5) Use Wi, Ui, Vi, and betai to calculate an improved estimate of b
    b = sum(Wi .* betai .* Vi) ./ sum(Wi .* betai .* Ui);

    % (6) Use the new b and repeat steps (3), (4), and (5) until successive
    % estimates of b agree within some desired tolerance tol
    b_diff = b - b_prev;
end

% (7) From this final value of b, together with the final x_bar and y_bar,
% calculate a from
a = y_bar - b .* x_bar;

% Goodness of fit
S = sum(Wi .* (yi - b.*xi - a).^2 );

