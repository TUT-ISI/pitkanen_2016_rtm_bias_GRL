function boxPlot_cell(D, varargin)

% boxPlot  Make a box and whisker plot
%
% Synopsis:  boxPlot_cell(D)
%            boxPlot_cell(D, x_group)
%            boxPlot_cell(D, x_group, ylimit)
%            boxPlot_cell(D, x_group, ylimit, width)
%            boxPlot_cell(D, x_group, ylimit, width, lineWidth)
%            boxPlot_cell(D, x_group, ylimit, width, lineWidth, verbose)
%            boxPlot_cell(D, x_group, ylimit, width, lineWidth, verbose, extra_info)
%
% Input:  D         = cell vector. D can contain NaNs which are ignored
%         x_group   = vector giving the box locations on x-axis
%         width     = width of each box.  Default: width=1
%         lineWidth = user-specfied line width for box and whiskers.  The
%                     median is drawn with twice this line thickness
%                     Default:  lineWidth = current default LineWidth
%                     lineWidth = get(gcf,'DefaultLineLineWidt');
%         verbose   = flag to turn on printing of summary data.
%                     Default:  verbose = false ==> no printing
%         extra_info  a cell variable. numbers or strings to write below
%                     each box
%
% Example 1:    Y{1} = randn(20,1); Y{2} = randn(20,1); Y{3} = randn(20,1);
%               boxPlot_cell(Y);
% Example 2:    Y{1} = randn(20,1); Y{2} = randn(20,1); Y{3} = randn(20,1);
%               boxPlot_cell(Y,[1 2 5]);
%
% Downloaded from MATLAB-Central on 2009-10-18, Original code by Shane Lin
% (no licence)
% Modified by Gerald Recktenwald, gerry@me.pdx.edu, on 2009-10-19
% Modified by Mikko Pitkänen, mikko.pitkanen@fmi.fi, June 2016
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 

% only want 7 inputs at most
numvarargs = length(varargin);
if numvarargs > 7
    error('boxPlot_cell requires at most 7 inputs');
end

%set the defaults for the optional arguments:
% {x_group, ylimit, width, lineWidth, verbose, extra_info}
optional_arg = {1:numel(D), [], 0.2, get(gcf,'DefaultLineLineWidth'), false, cell(size(D))};

% replace the ones defined by the user
optional_arg(1:numvarargs) = varargin;

% set easy-to-remember variable names
[x_group, ylimit, width, lineWidth, verbose, extra_info] = optional_arg{:};

% check inputs
if isempty(D);
   error('D cannot be empty in boxPlot_cell')
end

if numel(x_group) ~= numel(D);
   error('x_group and D must have the same number of elements in boxPlot_cell')
end

if isempty(width)
    width   = 0.2;
end

if isempty(lineWidth)   
    lineWidth = get(gcf,'DefaultLineLineWidth');
end

if isempty(verbose)
    verbose=false;
end

% calculate statistics and outlier for plotting
% Q = [highest_percentile; 75-perc; median; 25-perc; lowest-perc; number_observations]
[Q,Outliers] = getQuartiles(D);

ax = gca;
hold on;

%plot each bar according to statistics and add outliers
if 1==1
    for i = 1:numel(Q(1,:))
        Outliers{i} = [];      % optional: dont plot outliers
        drawOneBox(Q(:,i), Outliers{i}, width, lineWidth,x_group(i),ax);
        box on;
    end
    
% the same, but dont plot whiskers!
else
    
    for i = 1:numel(Q(1,:))
        Outliers{i} = [];      % optional: dont plot outliers
        data_wo_whiskers = [Q(2,i);Q(2,i);Q(3,i);Q(4,i);Q(4,i);Q(6,i);Q(7,i)];
        drawOneBox(data_wo_whiskers, Outliers{i}, width, lineWidth,x_group(i),ax);
        box on;
    end
end

%add data numbers and extra_info under boxes
if 1==1
    if ~(nargin<7 || isempty(extra_info))
        if isempty(ylimit)
            ylimit      = get(gca,'ylim'); 
        end
        yrange      = ylimit(2)-ylimit(1);
        ylimit(1)   = ylimit(1) - 0.1 * yrange;
        text_ypos   = ylimit(1)*ones(1,numel(D))+yrange*0.06;
        set(gca,'ylim',ylimit);
        text(x_group,text_ypos,num2str(Q(6,:)'),'HorizontalAlignment','center');


    end
end

%plot mean values
plot(x_group,Q(7,:),'r.');

%print verbose output
if verbose
    fprintf('\nQuartiles\n');
    fprintf('   j     Max       Q3       Q2       Q1      Min     NO\n');
    for j=1:size(Q,2)
        if isempty(Outliers{j})
            nout = 0;
        else
            nout = length(Outliers{j});
        end
        fprintf('%4d  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %4d\n',j,Q(:,j)',nout);
    end
end


%% ------------------------------------------------------------------------
function [Q,Outliers] = getQuartiles(D)
% getQuartiles  Find quartiles and outliers for a box-and-whisker plot
%
% Synopsis:  [Q,Outliers] = getQuartiles(D)
%
% Input: D  = cell vector
%            [m,n] = size(D) ==>  m values in each of n categories
%
% Output:  Q = matrix of quartiles.  size(Q) = [5,n].  Values in each
%              column of Q are in descending magnitude:
%                Q(1,j) = maximum value in category j
%                Q(2,j) = third quartile (q3)
%                Q(3,j) = second quartile (a.k.a. median, q2)
%                Q(4,j) = first quartile (q1)
%                Q(5,j) = minimum value in category j
%              The minimum (Q(5,j)) and maximum (Q(1,j)) are adjusted
%              to exclude outliers, see definitions below

Q = zeros(7,size(D,2));
Outliers = cell(size(D));

% -- Work through columns of D.  Code does not vectorize
%    because NaNs are removed, which causes the effective
%    length of each column to vary.
for j=1:size(D,2)
  x = D{j};
  x( isnan(x) ) = [];    %  remove NaNs
  
  if(isempty(x)); 
    Q(:,j) = [nan; nan; nan; nan; nan; nan; nan];
    continue;
  end;
  
  av = mean(x);

  % -- Compute Quartiles
  q2 = median(x);
  
  %medians:
  if(1==1)
    q1 = median( x(x<=q2) );   %  Q1 is median of lower half.
    q3 = median( x(x>=q2) );   %  Q3 is median of upper half.
  
  %or std's
  else
    x_std = std(x);  
    q1 = q2 - x_std;   %  Q1 is median of lower half.
    q3 = q2 + x_std;   %  Q3 is median of upper half.
  end
  IQR  = q3-q1;              %  Inter-quartile range
  qmax = max(x);            %  Tentative max
  qmin = min(x);            %  and min for boxplot

  % -- Outliers for boxplot are 1.5*IQR above Q3 and 1.5*IQR below Q1
  % this is manual switch so use this to select the desired definition for
  % the outliers
  if(1==1)
    klow  = x < (q1-1.5*IQR);
    khigh = x > (q3+1.5*IQR);
    
  %or you can use 10 and 90 percentiles for the whiskers
  elseif(1==1)
    x_sort = sort(x);
    lo_ind = floor(0.1*numel(x_sort));
    hi_ind  = ceil(0.9*numel(x_sort));
    if(lo_ind == 0); lo_ind = 1; end;
    klow  = x < x_sort(lo_ind);
    khigh = x > x_sort(hi_ind);
    
  %or you can use 5 and 95 percentiles for the whiskers
  elseif(1==1)
    x_sort = sort(x);
    lo_ind = floor(0.05*numel(x_sort));
    hi_ind  = ceil(0.95*numel(x_sort));
    if(lo_ind == 0); lo_ind = 1; end;
    klow  = x < x_sort(lo_ind);
    khigh = x > x_sort(hi_ind);
    
  %or 2*std's for the whiskers
  else
    x_std = std(x); 
    klow  = x < (q2 - 2*x_std);
    khigh = x > (q2 + 2*x_std);
  end
  
  iout = find( klow | khigh );
  if isempty(iout)
    Outliers{j} = [];
  else
    Outliers{j} = x(iout);
    if any(klow)
      qmin = min( x(~klow) );
    end
    if any(khigh)
      qmax = max( x(~khigh) );
    end
  end
  if(isempty(qmin) || isempty(qmax)); qmin = NaN; qmax=NaN; end;
  
  N = numel(x);
  Q(:,j) = [qmax; q3; q2; q1; qmin; N; av];
end


%% ------------------------------------------------------------------------
function drawOneBox(Q, Outliers, width, lineWidth, x, ax)

% drawBox  Create the box-and-whisker plot from quartile data
%
% Synopsis:  drawBox(Q, Outliers, lineWidth, width)
%
% Input:  Q = matrix of quartiles created by getQuartiles
%         Outliers = cell array of outliers for the box-and-whisker plot
%                    See getQuartiles
%         width = dimension used to specify the width of each box.  Use width=1
%         lineWidth = user-specfied line width for box and whiskers.  The
%                     median is drawn with twice this line thickness

% -- Set up a NaN-separated list of values for each box
%    19 points for each box & whisker pair.  20th point
%    is separatore between

%n = size(Q,2);   %  number of boxes.  Data for each box in column of Q
%n = 1;   %  number of boxes.  Data for each box in column of Q

%ib = (1:n);
ib = x;
%hwidth = (1-1/(1+n))/(1+9/(width+3));   %  half-width of box
hwidth = width;   %  half-width of box
qwidth = hwidth/2;                      %  quarter-width
ewidth = qwidth/2;                      %  eighth-width

%calculate notches
%       q2   - 1.57 .* (q3 - q1)     ./ sqrt(N)
no(1) = Q(3) - 1.57 .* (Q(2) - Q(4)) ./ sqrt(Q(6)); 
%       q2   + 1.57 .* (q3 - q1)         ./ sqrt(N)
no(2) = Q(3) + 1.57 .* (Q(2) - Q(4)) ./ sqrt(Q(6)); 

if(no(1) < Q(4));
    no(1) = Q(4);
end;

if(no(2) > Q(2));
    no(2) = Q(2);
end;

xb = zeros(23,1);           yb = xb;
xb(1)  = ib - ewidth;       yb(1) = Q(5);    %  min left
xb(2)  = ib + ewidth;       yb(2) = Q(5);    %  min right
xb(3)  = ib;                yb(3) = Q(5);    %  min middle
xb(4)  = ib;                yb(4) = Q(4);    %  Q1 middle
xb(5)  = ib + hwidth;       yb(5) = Q(4);    %  Q1 right
xb(6)  = ib + hwidth;       yb(6) = no(1);   %  Q1 notch right ---
xb(7)  = ib + qwidth;       yb(7) = Q(3);    %  Q2 right midwaist
xb(8)  = ib - qwidth;       yb(8) = Q(3);    %  Q2 left midwaist
xb(9)  = ib - hwidth;       yb(9) = no(2);   %  Q3 notch left  ---
xb(10)  = ib - hwidth;      yb(10) = Q(2);   %  Q3 left
xb(11)  = ib;               yb(11) = Q(2);   %  Q3 mid
xb(12) = ib;                yb(12) = Q(1);   %  max mid
xb(13) = ib - ewidth;       yb(13) = Q(1);   %  max left
xb(14) = ib + ewidth;       yb(14) = Q(1);   %  max right
xb(15) = ib;                yb(15) = Q(1);   %  max mid
xb(16) = ib;                yb(16) = Q(2);   %  Q3 mid
xb(17) = ib + hwidth;       yb(17) = Q(2);   %  Q3 right
xb(18) = ib + hwidth;       yb(18) = no(2);  %  Q3 notch right ---
xb(19) = ib + qwidth;       yb(19) = Q(3);   %  Q2 right midwaist
xb(20) = ib - qwidth;       yb(20) = Q(3);   %  Q2 left midwaist
xb(21) = ib - hwidth;       yb(21) = no(1);  %  Q1 notch left  ---
xb(22) = ib - hwidth;       yb(22) = Q(4);   %  Q1 left
xb(23) = ib;                yb(23) = Q(4);   %  Q1 mid

% -- Plot the basic box
plot(ax,xb,yb,'Linewidth',lineWidth);

% -- Add medians as thick lines
xm = zeros(2,1);  ym = xm;
xm(1) = ib - qwidth;   ym(1) = Q(3);
xm(2) = ib + qwidth;   ym(2) = Q(3);
hold('on');
plot(xm,ym,'Linewidth',2*lineWidth);

% -- Add outliers, if there are any
for j=1:size(Q,2)
  if any( Outliers )
    %nout = length(Outliers);
    plot(ax,ib,Outliers,'b.');
  end
end
hold('on');

