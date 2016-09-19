% A handy Matlab function that lets you print figures that all
% have approximately the same measures in centimeters independent of your
% computer and screen size, when you use the same resolution. Basically
% this function sets new values to default plotting and printing parameters
% and they can be easily reset by restarting Matlab or by running
% set_fig_properties again.
%
% Usage:
%     set_fig_properties(width,height,fsz,lw,msz)
%     width         figure width in [cm]
%     height        figure height in [cm]
%     fsz           font size
%     lw            line width
%     msz           markersize
%
% Example:
%     set_fig_properties(15,15,16,3,5)
%     figure
%     x = -pi:pi/10:pi;
%     y = tan(sin(x)) - sin(tan(x));
%     plot(x,y,'--rs')
%     print('-dpng','test1.png','-r300')
%     figure
%     plot(x,2*y,'--bs')
%     print('-dpng','test2.png','-r300')
%
% Run just once before plotting and all following figures will have
% the same properties. Warning: will close all open figures when executed.
%
% Code adopted from:
% http://dgleich.github.io/hq-matlab-figs/
% So credits go to:
% Tamara G. Kolda, Sandia National Laboratories and David F. Gleich,
% Purdue University
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

function [] = set_fig_properties(width,height,fsz,lw,msz)

% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all

% A guideline for different image media
%               Default     Paper       Presentation
% Width         5.6         varies      varies
% Height        4.2         varies      varies
% AxesLineWidth 0.5         0.75        1
% FontSize      10          8           14
% LineWidth     0.5         1.5         2
% MarkerSize	6           8           12

width   = width  / 2.54;    % cm -> in
height  = height / 2.54;    % cm -> in

% paper type is usletter by default, lets change to A4
set(0, 'DefaultFigurePaperType','A4')

% set the default line width to lw
set(0,'defaultLineLineWidth', lw);  

% set the default line marker size to msz
set(0,'defaultLineMarkerSize',msz); 

% set the default font size
set(0,'DefaultAxesFontSize',  fsz);

set(0,'DefaultAxesFontName', 'Arial')

% Set the default size for display
defpos = get(0,'defaultFigurePosition');

% set user defined height
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1) - width)/2;
bottom = (defsize(2) - height)/2;

% user defined height
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

close all
