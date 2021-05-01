function [frequency,amplitude,phase] = calcnaff(Y, Yp, varargin)
%CALCNAFF Computes NAFF decomposition for a phase space trajectory
% [nu amp phase] = calcnaff(Y, Yp, Win)
%
%  INPUTS
%  1. Y  - position vector
%  2. Yp -
%  3. WindowType  - Window type - 0 {Default} no windowing
%                                 1 Window of Hann
%                                 2 etc
%  4. nfreq - Maximum number of fundamental frequencies to search for
%             10 {Default}
%  5. debug - 1 means debug flag turned on
%             0 {Default}
%
%  Optional Flags
%  'Debug' - turn on deubbing flag
%  'Display' - print ou results
%  'Hanning' - select Window of Hann, WindowType = 1
%  'Raw' or 'NoWindow' - select Window of Hann, WindowType = 0
%
%  OUTPUTS
%  1. frequency - frequency vector with sorted amplitudes
%                 by default the algorithm try to compute the 10 first fundamental
%                 frequencies of the system.
%  2. amplitude - amplitudes associated with fundamental frequencies
%  3. phase - phases associated with fundamental frequencies
%
%  NOTES
%  1. Mimimum number of turns is 64 (66)
%  2. Number of turn has to be a multiple of 6 for internal optimization
%  reason and just above a power of 2. Example 1026 is divived by 6 and
%  above 1024 = pow2(10)
%
%  Examples
%  NT = 9996; % divided by 6
%  simple quasiperiodic (even period) motion 
%  y =2+0.1*cos(pi*(0:NT-1))+0.00125*cos(pi/3*(0:NT-1));
%  yp=2+0.1*sin(pi*(0:NT-1))+0.00125*sin(pi/3*(0:NT-1));
% 
%  [nu ampl phase] = calcnaff(y,yp,1); % with windowing


% Written by Laurent S. Nadolski
% April 6th, 2007
% Modification September 2009: 
%  test if constant data or nan data

% BUG in nafflib: returns nan even if valid data. Number of try
nitermax = 10;

% Flag factory
[wraw1,args]=getflag(varargin,'Raw'); %#ok<ASGLU>
[wraw2,args]=getflag(args,'NoWindow'); %#ok<ASGLU>
[whann,args]=getflag(args,'Hanning');
[dbg,args]=getflag(args,'Debug');
[DisplayFlag,args]=getflag(args,'Display');
[WindowType,nfreq,DebugFlag]=getargs(args,0,10,double(dbg));
if whann, WindowType=1; end


% Test wether nan or constant data
if any(isnan(Y(1,:)))
    fprintf('Warning Y contains NaNs\n');
    frequency =NaN; amplitude = NaN;  phase = NaN;
elseif any(isnan(Y(1,:)))
    fprintf('Warning Yp contains NaNs\n');
    frequency =NaN; amplitude = NaN;  phase = NaN;
elseif (mean(Y) == Y(1) && mean(Yp) == Yp(1))
    fprintf('Warning data are constant\n');
    frequency = 0; amplitude = 0;  phase = 0;
else % Frequency map analysis
    [frequency,amplitude,phase] = nafflib(Y, Yp, WindowType,nfreq,DebugFlag);
    %It seems there is a bug in nafflib, something returns nan even for valid data 
    niter = 0;
    while any(isnan(frequency)) && (niter < nitermax)
        pause(2);
        fprintf('Warning Nan returned by NAFF (x%d)\n', niter);
        niter = niter +1;
        [frequency,amplitude,phase] = nafflib(Y, Yp, WindowType,nfreq,1); % add debugging
    end
        
    if DisplayFlag
        fprintf('*** Naff run on %s\n', datestr(clock))
        for k = 1:length(frequency)
            fprintf('nu(%d) =% .15f amplitude = % .15f phase = % .15f \n', ...
                k,frequency(k), amplitude(k),phase(k));
        end
    end
end
