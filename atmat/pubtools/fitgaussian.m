function varargout = fitgaussian(varargin)

% GAUSSIAN_PARAM FITERR GAUSSFIT SIGERROR]= FITGAUSSIAN(DATA,[property_value_pair]);
% 
% DATA is a 1D vector to which we want to fit a gaussian profile. The
% function will return a structure with various fitted parameters.
% SCALEFACTOR is optional and is applied in the horizontal axis.
%
% Property value pairs:
%
%  scale    : 1 (default)
%  plot     : 0 no plots (default), 1 plots fits against data
%
% Initial Fit parameters are fit automatically unless specified.
%  integral : estimated sum integral of the function
%  mean     : estimate of the mean of the gaussian 
%  sigma    : estimated sigma of the gaussian
%  DC       : known DC component to remove (set to zero to not fit)
%  bg_grad  : known background gradient to fit (set to zero to not fit)
%  Assym    : known gaussian assymetry (set to zero to not fit)
%
% FITERR    Final fit error returned by the internal error fuction 
% GAUSSFIT  Final fit function
% SIGERROR  Error estimate of the sigma fit in relative terms. Multiply by
%           100 to get in percent.
%
% Original script by R. Dowd
% In functional format by E. Tan 31/07/2009

DEBUG = 0;
ASYM_FLAG = 0;

[reg, prop] = parseparams(varargin);

if nargin > 0
    data = reg{1};
else
    disp('No data to work with');
    return
end
data = data(:)';

% defaults
scalefac = 1;
plotfig = 0;
autofitDC = 1;
autofitGrad = 1;
autofitAssym = 1;
userx = [];
for i=1:length(prop)/2
    ind = (i-1)*2+1;
    switch lower(prop{ind})
        case 'scale'
            scalefac = prop{ind+1};
        case 'plot'
            plotfig = prop{ind+1};
        case 'integral'
            guess_area = prop{ind+1};
        case 'mean'
            guess_mu = prop{ind+1};
        case 'sigma'
            guess_sigma = prop{ind+1};
        case 'dc'
            DCfit = prop{ind+1};
            autofitDC = 0;
        case 'bg_grad'
            linefit = prop{ind+1};
            autofitGrad = 0;
        case 'assym'
            Assym = prop{ind+1};
            autofitAssym = 0;
        case 'x'
            userx = prop{ind+1};
    end
end

% percentage of at the start and end of the data set assumed to be
% representative of the background noise.
pcbackground = 0.1;

datasize = length(data);

if isempty(userx)
    x=[1:datasize]*scalefac;
else
    x = userx;
end

% Fit a gaussian; first find some starting parameters

% Used to guess the DC component. Assume flat for the first 10% of data
% points.
startDC = mean(data(1:fix(datasize*pcbackground)));
endDC =  mean(data(end-fix(datasize*pcbackground):end));


if ~exist('guess_area','var')
    % Guess area AUTOMATICALLY
%     guess_area = sum(data) - 0.5*(startDC+endDC)*datasize;
%     guess_area = guess_area*scalefac;
    guess_area = sum((data(2:end) + data(1:end-1)).*diff(x)/2) - 0.5*(startDC+endDC)*(x(end)-x(1));
    guess_area = guess_area;
end

if ~exist('guess_mu','var')
    % Guess the center of mass AUTOMATICALLY
    [~, maxind] = max(data);
    guess_mu = x(maxind);
end

if ~exist('guess_sigma','var')
    % Guess sigma in pixels AUTOMATICALLY
    maxval = max(data);
    indices = find(data > (maxval+((startDC+endDC)/2) )/2);
    guess_sigma = (x(indices(end)) - x(indices(1)))/2.3;
    guess_sigma = guess_sigma;
end

% So far everything has been calculated in units of data points. Apply
% scaling factor here.
fixedvals = [NaN NaN NaN];

Starting(1) = guess_area;
Starting(2) = guess_mu;
Starting(3) = guess_sigma;
if autofitDC
    Starting(end+1) = startDC;
else
    fixedvals(1) = DCfit;
end
if autofitGrad
    % Guess if there is a background gradient. Again assuming first and last
    % 10% of data set is "background".
    Starting(end+1) = -(endDC-startDC)/datasize*scalefac;
else
    fixedvals(2) = linefit;
end
if autofitAssym
    % Initial Assymetry factor
    Starting(end+1) = 0;
else
    fixedvals(3) = Assym;
end


if DEBUG
    options = optimset('Display','iter','MaxIter',1500,'TolX',1e-6,'TolFun',1e-10);
else        
    options = optimset('Display','off','MaxIter',1500,'TolX',1e-6,'TolFun',1e-10);
end

[Estimates fval] = fminsearch(@myfit,Starting,options,x,data,fixedvals);


fitparam.xdata = x;
fitparam.rawdata = data;
fitparam.area = Estimates(1);
fitparam.mu = Estimates(2);
fitparam.sigma = Estimates(3);
i = 1;
if autofitDC
    fitparam.DC = Estimates(3+i);
    i = i + 1;
else
    fitparam.DC = fixedvals(1);
end
if autofitGrad
    fitparam.bg_gradient = Estimates(3+i);
    i = i + 1;
else
    fitparam.bg_gradient = fixedvals(2);
end
if autofitAssym
    fitparam.Assym_factor = Estimates(3+i);
else
    fitparam.Assym_factor = fixedvals(3);
end
fitparam.final_fit_val = fval;



varargout{1} = fitparam;
if nargout > 1
    varargout{2} = fval;
end
if nargout > 2
    gaussianfit = ones(size(data));
    for i = 1:datasize
        c = x(i);
        gaussianfit(i) = fitparam.area * exp(-0.5*((c-fitparam.mu)./((1+sign(c-fitparam.mu)*fitparam.Assym_factor)*fitparam.sigma)).^2) / sqrt(2*pi*fitparam.sigma^2) + fitparam.bg_gradient*c + fitparam.DC;
    end
    varargout{3} = gaussianfit;
end
if nargout > 3
    % Calculate error in sigma
    er=[];
    errscale = ones(size(Estimates));
    for perc=0.90:0.001:1.1;
        errscale(3) = perc; % change the sigma value and see the fit.
        er(end+1) = myfit(Estimates.*errscale,x,data,fixedvals);
    end
    % Normalise
    er = er./min(er);
    % threshold a 5% change in the error function;
    ind = find(er<1.05);
    perc = 0.90:0.001:1.1;
    % percerror
    sigmaerror = (perc(ind(end))-1);
    
    varargout{4} = sigmaerror;
end

if DEBUG || plotfig
    gaussianfit = ones(size(data));
    for i = 1:datasize
        c = x(i);
        gaussianfit(i) = fitparam.area * exp(-0.5*((c-fitparam.mu)./((1+sign(c-fitparam.mu)*fitparam.Assym_factor)*fitparam.sigma)).^2) / sqrt(2*pi*fitparam.sigma^2) + fitparam.bg_gradient*c + fitparam.DC;
    end
    figure(233);
    plot(x,data, '.-r');
    hold on;
    plot(x,gaussianfit, '-b');
    hold off;
%     title(sprintf('Fitting STD error %g (\\sigma = %f)',...
%         fittingerror(end),fitsigmas(end)));
end




function sse=myfit(params, x, Dist, fixedvals)

% if length(params) > 0
    Afit = params(1);
% else
%     Afit = 1;
% end
% if length(params) > 1
    mufit = params(2);
% else
%     mufit = length(x)/2;
% end
% if length(params) > 2
    sigmafit = params(3);
% else
%     sigmafit = length(x)/10;
% end

i = 1;
if isnan(fixedvals(1))
    DCfit = params(3+i);
    i = i + 1;
else
    DCfit = fixedvals(1);
end

if isnan(fixedvals(2))
    linefit = params(3+i);
    i = i + 1;
else
    linefit = fixedvals(2);
end

if isnan(fixedvals(3))
    Asym = params(3+i);
else
    Asym = fixedvals(3);
end

fittedcurve = Afit * exp(-0.5*((x-mufit)./((1+sign(x-mufit)*Asym)*sigmafit)).^2) / sqrt(2*pi*sigmafit^2) + (linefit*x) + DCfit;

sse = sum((fittedcurve - Dist).^2)/length(Dist);