function violin( x, y, varargin)
%% VIOLIN Violin plot                               
% Input
%   x           single numeric value
%   y           numeric vector
%
%   Optional input (key-value pairs, case-insensitive):
%   'n'         Scalar, number of values of non-parametric kernel
%               Default = 1000
%   'facecolor' Color, according to matlab plot color properties
%               Default = [0 0 0] = black
%   'facealpha' Transparancy, where 0=fully transparent, 1=fully opaque.
%               Default = 0.5
%   'withmdn'   Logical, if true median is plotted as a line. 
%               Default = false
%   'linestyle' Linestyle, according to matlab plot linestyle properties
%               Default = '-'
%   'linewidth' Linewidth, according to matlab plot linestyle properties
%               Default = 1
%   'linecolor' Linecolor, according to matlab plot linestyle properties
%               Default = 'k'
%   'style'     1 = opaque outline + tranparent violin patch only
%               2 = transparent violin patch only
%               3 = opaque outline only
%               Default = 1
%   'side'      'both'  = plot two-sided violin (Default)
%               'left'  = plot only left-sided violin
%               'right' = plot only right-sided violin
%   'rotation'  orientation of violin plots
%               'vertical'   = violin stands upright
%               'horizontal' = violin lays down
%   'scaling'   scalar, scaling of width of violin, in x coordinates.
%               numerical value
%               Default = 1
%   'cutoff'    cutoff density
%               Default = 1e-3
%   'kernel'    The type of kernel smoother to use, chosen from 
%               'normal'  'box', 'triangle', and 'epanechnikov'.
%               Default = 'normal'
%   'kernelwidth' scalar, for the bandwidth of the kernel smoothing window.
%                 The default is optimal for estimating normal densities,
%               but you may want to choose a smaller value to reveal
%               features such as multiple modes.
%
% EXAMPLE
% figure(1); clf;
% rng(17371218);
% violin(1,randn(1,1000), ...
%        'KernelWidth', 0.5, ...
%        'withmdn', 1);
% violin(0,randn(1,1000) + 6, ...
%        'KernelWidth', 0.5, ...
%        'Rotation',    'horizontal', ... 
%        'side',        'right', ...
%        'facecolor',   'r', ...
%        'cutoff',      1e-3, ...
%        'scaling',     4);
% 
%                                                               J.H.F. 2019
%-------------------------------------------------------------------------%
%% Parse variable input arguments                   
for k = 1:2:length(varargin)
    switch lower(varargin{k})
        case 'n'
            n = varargin{k+1};
        case 'facecolor'
            facecolor = varargin{k+1};
        case 'facealpha'
            facealpha = varargin{k+1};
        case 'withmdn'
            withmdn = varargin{k+1};
        case 'linestyle'
            linestyle = varargin{k+1};
        case 'linewidth'
            linewidth = varargin{k+1};
        case 'linecolor'
            linecolor = varargin{k+1};
        case 'style'
            style = varargin{k+1};
        case 'side'
            side = lower(varargin{k+1});
            if ~ismember(side,{'both','left','right'})
                error('%s is not a valid value for ''side'' argument in violin',varargin{k});
            end
        case 'rotation'
            rotation = lower(varargin{k+1});
            if ~ismember(rotation,{'vertical','horizontal'})
                error('%s is not a valid value for ''rotation'' argument in violin',varargin{k});
            end
        case 'scaling'
            scaling = varargin{k+1};
        case 'cutoff'
            cutoff = varargin{k+1};
        case 'kernel'
            kernel = lower(varargin{k+1});
            if ~ismember(kernel,{'normal','box','triangle','epanechnikov'})
                error('%s is not a valid value for ''kernel'' argument in violin',varargin{k});
            end
        case 'kernelwidth'
            kernelwidth = varargin{k+1};
        otherwise
            error('''%s'' is not a valid key input argument for violin',varargin{k});
    end
end
%% Set default parameter values                     
if ~exist('n','var') || isempty(n)
    n = 1000;
end
if ~exist('facecolor','var') || isempty(facecolor)
    facecolor = 'k';
elseif isnumeric(facecolor) && length(facecolor)==1
    facecolor = ones(1,3)*facecolor;
end
if ~exist('facealpha','var') || isempty(facealpha)
    facealpha = 0.5;
end
if ~exist('withmdn','var') || isempty(withmdn)
    withmdn = false;
end
if ~exist('linestyle','var') || isempty(linestyle)
    linestyle = '-';
end
if ~exist('linewidth','var') || isempty(linewidth)
    linewidth = 1;
end
if ~exist('linecolor','var') || isempty(linecolor)
    linecolor = 'k';
elseif isnumeric(linecolor) && length(linecolor)==1
    linecolor = ones(1,3)*linecolor;
end
if ~exist('style','var') || isempty(style)
    style = 1;
end
if ~exist('side','var') || isempty(side)
    side = 'both';
end
if ~exist('rotation','var') || isempty(rotation)
    rotation = 'vertical';
end
if ~exist('scaling','var') || isempty(scaling)
    scaling = 1;
end
if ~exist('cutoff','var') || isempty(cutoff)
    cutoff = 1e-3;
end
if ~exist('kernel','var') || isempty(kernel)
    kernel = 'normal';
end
if ~exist('kernelwidth','var') || isempty(kernelwidth)
    kernelwidth = [];
end
switch style
    case 1
        % do nothing everything is set with current parameter settings
    case 2
        linestyle = 'none';
    case 3
        facealpha = 0;
    otherwise
        error('Style can only take numbers 1, 2 or 3');
end
%% Fit y values (nonparametric kernel-smoothing)    
% fit distribution
pd = fitdist(y(:),'kernel','kernel',kernel,'width',kernelwidth);
% ypdf
ymin = min(y)-(max(y)-min(y));
ymax = max(y)+(max(y)-min(y));
ypdf_tmp = linspace( ymin, ymax, n );
% xpdf
xpdf_tmp  = pdf(pd,ypdf_tmp);
if ischar(scaling) && strcmp(scaling,'off')
    scaling = max(xpdf_tmp);
end
% adjust amplitude to max width
xpdf_norm = xpdf_tmp * scaling;
% length of median line
if withmdn
    mdnlength = xpdf_norm(find(ypdf_tmp>nanmedian(y),1,'first')-1);
end
%% Trim at cutoff                                   
include = xpdf_tmp > cutoff;
xpdf_inc = xpdf_norm(include);
ypdf_inc = ypdf_tmp(include);
%% Plotting values                                  
xpdfR = x + xpdf_inc;
xpdfL = x - xpdf_inc;
switch side
    case 'both'
        xpdf = [xpdfR(:)' fliplr(xpdfL(:)')];
        ypdf = [ypdf_inc(:)' fliplr(ypdf_inc(:)')];
    case 'left'
        xpdf = fliplr(xpdfL(:)');
        ypdf = fliplr(ypdf_inc(:)');
        xpdf = [x xpdf(:)' x];
        ypdf = [ypdf(1) ypdf(:)' ypdf(end)];
    case 'right'
        xpdf = [x xpdfR(:)' x];
        ypdf = [ypdf_inc(1) ypdf_inc(:)' ypdf_inc(end)];
end
switch rotation
    case 'vertical'
        % do nothing everything is already set
    case 'horizontal'
        xpdfVERT = xpdf;
        ypdfVERT = ypdf;
        ypdf = xpdfVERT;
        xpdf = ypdfVERT;
end
%% Plot                                             
violinplot = patch( xpdf, ypdf, facecolor);
% beautify
set( violinplot, 'EdgeColor', linecolor );
set( violinplot, 'LineStyle', linestyle );
set( violinplot, 'LineWidth', linewidth );
set( violinplot, 'FaceAlpha', facealpha );
% plot median
if withmdn
    hold on
    switch rotation
        case 'vertical'
            switch side
                case 'both'
                    medianplot = plot( x+[-1 1].*mdnlength, ones(1,2)*nanmedian(y) );
                case 'left'
                    medianplot = plot( x+[-1 0].*mdnlength, ones(1,2)*nanmedian(y) );
                case 'right'
                    medianplot = plot( x+[ 0 1].*mdnlength, ones(1,2)*nanmedian(y) );
            end
        case 'horizontal'
            switch side
                case 'both'
                    medianplot = plot( ones(1,2)*nanmedian(y), x+[-1 1].*mdnlength );
                case 'left'
                    medianplot = plot( ones(1,2)*nanmedian(y), x+[-1 0].*mdnlength );
                case 'right'
                    medianplot = plot( ones(1,2)*nanmedian(y), x+[ 0 1].*mdnlength );
            end
    end
    
    % beautify
    set(medianplot, 'Color', linecolor );
    set(medianplot, 'LineWidth', linewidth );
    
end
