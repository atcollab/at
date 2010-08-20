function varargout = freqsearch3(data, varargin)
% =========================================================================
% Find the frequency terms in the data set with the use of filters and FFT
% or a "search" algorithm (slower but more accurate). The function returns
% the number of oscillations per unit time where the unit time is define by
% DT. Eg If DT is in seconds, then freq is the number of oscillations per
% second. If DT is the number of turns, then freq is the number of
% oscillations per turn.
%
% [freq(s) amplitude(s) eigenvector(s) time_vec_used] = ...
%          FREQSEARCH(DATA, [FILTER, METHOD, DT, ORDER, RANGE, TOLERANCE, windowfraction])
%
% DATA   : input data, can be complex.
% FILTER : 'hanning','none' (default)
% METHOD : 'fft' (default),'search','spectrum'
% DT     : timestep between each data point. (default: 1)
%
% The options below are only applicable to the 'search' method.
%
% ORDER (vector) : number of frequency terms to extract. Ordered by relative
%                  strength. So (default: [1])
% RANGE          : frequency range over which to scan. (default: [0 Inf])
% TOLERANCE      : Search until freq_(n) - freq_(n-1) < TOLERANCE (default:
%                  1e-10)
% windowfraction : How wide should the subsequent search range should be.
%                  (default: 0.03)
%
% Examples:
% >> [f a] = freqsearch(data,'hanning','search',1,[1 2 3])
%
% 24/01/2006
% Eugene
% v2.0 - fft and "search" method combines to increase the speed in which
%        one can analyse the freqency components.
%      - Use of filters, modified hanning type filter only at the moment.
%      - Multiple orders for comparisons of multiple resonant frequencies.
%        The frequency is ordered by amplitude/strength. The first
%        frequency being the dominant one followed by the second strongest.
% 19/08/2010 Eugene: added 'spectrum' option to return a spectrogram. When
%                    in this mode, frequency and amplitude will be vectors
%                    that represent the spectrogram. Range will need to be
%                    specified and has to be a vector at the frequencies of
%                    interest.

DEBUG = false;

%======================================================================
% Parse input
% Set some defaults (mainly for the search method)
% Min and Max number of iterations for search method
MAXIT = 23;

% What filter to use
Nparam = 1;
if nargin >= Nparam + 1 & ischar(varargin{Nparam})
    switch lower(varargin{Nparam})
        case 'hanning'
            filter = 'hanning';
        case 'none'
            filter = 'none';
        otherwise
            error(sprintf('Unknown filter option %s',varargin{Nparam}));
    end
else
    filter = 'none';
end

% What method to use
Nparam = Nparam + 1;
if nargin >= Nparam + 1 & ischar(varargin{Nparam})
    method = lower(varargin{Nparam});
else
    method = 'fft';
end  

% Time step between each sample
Nparam = Nparam + 1;  
if nargin >= Nparam + 1
    dt = varargin{Nparam};
else
    dt = 1;
end

% Number of terms to extract
Nparam = Nparam + 1;
if nargin >= Nparam + 1
    order = varargin{Nparam};
else
    order = 1;
end  

% Tune range
Nparam = Nparam + 1;  
if nargin >= Nparam + 1
    range = varargin{Nparam};
else
    range = [0 0.5];
end

Nparam = Nparam + 1;  
if nargin >= Nparam + 1
    tolerance = varargin{Nparam};
else
    tolerance = 1e-10;
end

% Determines how much to zoom in when narrowing the search.
% Depending on the tolerance, the optimal value for the windowfraction
% changes. However 4% seems good enough.
Nparam = Nparam + 1;  
if nargin >= Nparam + 1
    windowfraction = varargin{Nparam};
else
    windowfraction = 0.04;
end

if DEBUG
    fprintf('Options selected: filter(%s) method(%s) dt(%11.3e)\n',...
        filter, method, dt);
    fprintf('                  order(%d) range(%f %f) tolerance(%11.3e) windowfraction(%f)\n',...
        order(end), range(1), range(2), tolerance, windowfraction);
end
% Finshed parsing input
%======================================================================

% Define variables
% Define the time or running parameter against which to calculate the
% freqrange.
neval = length(data);
T2 = dt*(neval-1)/2;  %-T/2 --> T/2
t = [-T2:dt:T2]';

eigenvec = zeros(neval,max(order));
orthvec = zeros(neval,max(order));
orthvec_ = zeros(neval,max(order));
a = zeros(1,max(order));
nu = zeros(1,max(order));

% Ensure that data and t are column vectors;
data = reshape(data,neval,1);
datareal = isreal(data);
%======================================================================

% Remove any DC component in the signal
% data = data - 0.5/T2*local_midpointsum(data);

% What filter to apply to the data
if DEBUG; disp('Calculating filter'); end;
usefilter = 0;
switch filter
    case 'hanning'
        % Window function that increases hight of the fundamental peak to make it
        % easier to pickout.
        p = 1; % cosine window order
        kai = complex( 2^p*(factorial(p))^2*(1+cos(pi*(t/T2))).^p/factorial(2*p) );
        usefilter = 1;
end
% Finished applying filter
%======================================================================

% What method to use
if DEBUG; disp('Starting calculation'); end;
switch method
    case 'fft'
        [nu a] = local_calculate_with_fft(data,dt,range);
        order = 1;

    case 'search'

        for k=1:max(order)
            % Start the frequency search using a two step approach, first
            % use the FFT to get a coarse measurement of the frequency
            % followed by the correlation analysis to get a more accurate
            % measure of the dominant frequency component.

            % FFT
            if usefilter
                prelim_freq = local_calculate_with_fft(data.*kai,dt,range);
            else
                prelim_freq = local_calculate_with_fft(data,dt,range);
            end

            % Will scan this range of frequencies.
            freqrange = local_find_new_range(prelim_freq,range(2),range(1),windowfraction);

            % Some initial variables. Start with some guess at the
            % frequency, mainly for the first difference calculation.
            % This is the power spectrum/frequency scan.
            psi = zeros(1,length(freqrange));
            freq = zeros(1,MAXIT);

            omega_prev = median(freqrange);
            difference = 1;

            for j=1:MAXIT
                % Do the integral that calculates the average <f(t), e^i*freqrange*t>. Not
                % including multiplication of some factors like dt, since we only
                % need to find where psi is a maximum and extract the corresponding
                % freqrange. Vectorising this loop does not help,
                % evaluated already.
                if usefilter
                    psi = local_psi_integral(data.*kai,t,freqrange);
                else
                    psi = local_psi_integral(data,t,freqrange);
                end

                if j >= 1 && j <=1 && DEBUG
                    figure; plot(freqrange,abs(psi));
                    xlabel('freq / Frequency'); ylabel('Arb. Units');
                end

                % Calculate the value of freqrange for the maximum psi.
                [maxpsi maxind] = max(psi(:));
                freq(j) = freqrange(maxind);

                difference = abs(freq(j) - omega_prev);
                if difference < tolerance
                    if DEBUG; fprintf('Difference less than specified tolerance. j=%d\n',j); end;
                    break;
                else
                    omega_prev = freq(j);
                end

                % Find new range to seach, zoom in.
                freqrange = local_find_new_range(freq(j),freqrange(end),freqrange(1),windowfraction);

                psi = zeros(size(freqrange));
            end
            if DEBUG; fprintf('FREQ = %20.10e\n',freq(1:j)); end;
            % Orthogonal projection to determine the coeffients. Since
            % e^i*2pi*freq*t.

            eigenvec(:,k) = exp(complex(0,2*pi*freq(j).*(t)));
            % Orthogonalize
            %             sumprojections = zeros(neval,1);
            %             for ii=1:k-1
            %                 sumprojections = sumprojections + dot(eigenvec(:,k),orthvec(:,ii))/dot(orthvec(:,ii),orthvec(:,ii))*orthvec(:,ii);
            %             end
            %             orthvec(:,k) = eigenvec(:,k) - sumprojections;

            a(k) = ((0.5/T2)*local_midpointsum(data.*conj(eigenvec(:,k))))*dt;
            %             a(k) = (0.5/T2)*maxpsi;
            nu(k) = freq(j);

            % Subtract the component from 'f' function.
            data = data - a(k)*eigenvec(:,k);
        end
    case 'spectrum'
        % Return the power spectrum
        if usefilter
            psi = local_psi_integral(data.*kai,t,range);
        else
            psi = local_psi_integral(data,t,range);
        end
        if datareal
            nu = psi*dt/T2;
        else
            nu = 0.5*psi*dt/T2;
        end
        order = 1:length(psi);
                
    otherwise
        error(sprintf('Unknown method option %s',varargin{Nparam}));
end

varargout{1} = nu(order);
if nargout > 1
    if datareal
        % With only real data the returned amplitudes should also be real.
        % And the factor 2 is needed here but not quite sure why just yet.
        varargout{2} = 2*abs(a(order));
    else
        varargout{2} = a(order);
    end
end
if nargout > 2
    varargout{3} = eigenvec(:,order);
end
if nargout > 3
    varargout{4} = t;
end

% temp = freq(find(freq ~= 0));
% varargout{1} = temp(end)/(2*pi);
% DEBUG
% fprintf('%i    %17.15g\n',j, difference);


function fctnsum = local_midpointsum(fctn)
% Vectorise the midpoint method of integrating a numerical function.
% f_n      = [fctn 0];
% f_nplus1 = [0 fctn];  % shift all the numbers one "space" to the right.
% midpoints = 0.5*(f_n + f_nplus1);

midpoints = 0.5.*(fctn(1:end-1) + fctn(2:end));
fctnsum = sum(midpoints);


function psi = local_psi_integral(data,t,freqrange)


midpoints = zeros(1,length(t)-1);
omegarange = -freqrange*2*pi;
fctn = zeros(1,length(t));
psi = zeros(1,length(omegarange));

for k=1:length(freqrange)
    fctn = data.*exp(complex(0,omegarange(k)).*t);
%     fctn = data.*complex(cos(freqrange(k).*t),-sin(freqrange(k).*t));
    midpoints = 0.5.*(fctn(1:end-1) + fctn(2:end));
    psi(k) = sum(midpoints);
end


function [freq varargout] = local_calculate_with_fft(data,dt,range)
% Find peak with FFT within "range" of frequencies.
% Auto calculate number of points to calculate fft. Use maximum
nn = [4:15];
ind = max(find(2.^nn - length(data) < 0));
Nfft = 2^nn(ind);

% Calculate FFT and the power spectrum
yy = fft(data,Nfft);
Pyy = yy.*conj(yy);

% Corresponding frequency range
f = 1/(dt*Nfft).*(0:Nfft/2);
ii = find(f > range(1) & f < range(2));

% Find peak
[maxval maxind] = max(Pyy(ii));
freq = f(ii(maxind));
if nargout > 1
    varargout{1} = abs(yy(ii(maxind)))/(Nfft/2);
end


function freqrange = local_find_new_range(centre,upper,lower,windowfraction)

% Find new range to seach, zoom in.
new_width = (upper - lower)*windowfraction;
% minimum frequency separation
min_freq_sep = new_width/500;

if centre == lower
    lowerbound = centre - new_width*2/windowfraction;
    upperbound = centre;
    
    freqrange = [lowerbound:(upperbound-lowerbound)/100:upperbound];
elseif centre == upper
    lowerbound = centre;
    upperbound = centre + new_width*2/windowfraction;
    
    freqrange = [lowerbound:(upperbound-lowerbound)/100:upperbound];
else
    lowerbound = centre - new_width;
    upperbound = centre + new_width;
    
%     num = 15;
%     scalefactor = (2*new_width - min_freq_sep*(num+1))/(num+2);
%     freqrange = lowerbound + cumsum((1 + cos(0:2*pi/num:2*pi))*scalefactor + min_freq_sep);
    
    scalefactor = (2*new_width - min_freq_sep*16)/17;
    freqrange = lowerbound + cumsum((1 + cos(0:2*pi/15:2*pi))*scalefactor + min_freq_sep);
end

% freqrange = freqrange*2*pi;


