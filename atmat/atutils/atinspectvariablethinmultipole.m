function [listpola,listpolb] = atinspectvariablethinmultipole(element,varargin)
%        [pola_array,polb_array] = atinspectvariablethinmultipole(element)
%
%        Get the atvariablethinmultipole polynom values per turn.
%
%        Keyword arguments
%            turns(int): Default 1. Number of turns to calculate.
%            T0(float): revolution time in seconds. Use only in SINE mode.
%            tparticle(float): Default 0. Time of the particle in seconds.
%
%        Returns
%            Two cell arrays with the values of PolynomA and PolynomB per turn.
    p = inputParser;
    addOptional(p,'turns',1,@(x) isnumeric(x) && isscalar(x) && (x > 0));
    addOptional(p,'T0',0,@(x) isnumeric(x) && isscalar(x) && (x > 0));
    addOptional(p,'tparticle',0,@(x) isnumeric(x) && isscalar(x) && (x > 0));
    parse(p,varargin{:});
    par = p.Results;
    
    turns = par.turns;
    T0 = par.T0;
    tparticle = par.tparticle;
    mode = element.Mode;
    
    switch mode
      case 0
        % revolution time
        timeoffset = T0 + tparticle;
      case 2
        % particle time
        timeoffset = tparticle;
      otherwise
        timeoffset = 0;
    end
    
    if isfield(element,'Ramps')
      ramps = element.Ramps;
    else
      ramps = 0;
    end
    
    if isfield(element,'Periodic')
      periodic = element.Periodic;
    else
      periodic = 0;
    end
    
    maxorder = element.MaxOrder;
    
    pola = nan(1,maxorder+1);
    polb = nan(1,maxorder+1);
    
    listpola = cell(1, turns);
    listpolb = cell(1, turns);
    
    for turn = 0:turns-1
      for order = 0:maxorder
        if isfield(element,'AmplitudeA')
          pola(order+1) = get_pol(element, 'A', ramps, mode,  ...
              timeoffset * turn, turn, order, periodic);
        end
        if isfield(element,'AmplitudeB')
          pola(order+1) = get_pol(element, 'B', ramps, mode, ...
              timeoffset * turn, turn, order, periodic);
        end
      listpola{turn+1} = pola;
      listpolb{turn+1} = polb;
      end
    end

end

function ampt = get_amp(amp, ramps, t)
%    get_amp(ele, amp, ramps, t)
%
%    get_amp returns the input value `amp` when ramps is False.
%    If ramps is True, it returns a value linearly interpolated
%    accoding to the ramping turn.
    ampt = amp;
    if length(ramps) == 4
        if t <= ramps(0)
            ampt = 0.0;
        elseif t <= ramps(1)
            ampt = amp * (t - ramps(0)) / (ramps(1) - ramps(0));
        elseif t <= ramps(2)
            ampt = amp;
        elseif t <= ramps(3)
            ampt = amp - amp * (t - ramps(2))/ (ramps(3) - ramps(2));
        else
            ampt = 0.0;
        end
    end
end


function ampout = get_pol(element , ab, ramps, mode, t, turn, order, periodic)
%   get_pol(element , ab, ramps, mode, t, turn, order, periodic)
%
%   get_pol returns the polynom a or b for a given mode, turn, order,
%   time and periodicity.
    allamp = element.(strcat("Amplitude",ab) );
    amp = allamp(order+1);
    ampout = 0;

    if amp == 0
        return
    end

    % get the ramp value
    ampoutaux = get_amp(amp, ramps, turn);

    switch mode 
      case 0
        % sin mode parameters
        whole_sin_above = element.(strcat("Sin",ab,"above"));
        freq = element.(strcat("Frequency",ab));
        ph = element.(strcat("Phase",ab));
        sinval = sin(2 * pi * freq * t + ph);
        if sinval >= whole_sin_above
            ampout = ampoutaux * sinval;
        else
            ampout = 0;
        end
      case 1
        ampout = nan;
      case 2
        nsamples = element.(strcat("NSamples",ab));
        if periodic || turn < nsamples
            func = element.(strcat("Func",ab));
            funcderiv1 = element.(strcat("Func",ab,"deriv1") );
            funcderiv2 = element.(strcat("Func",ab,"deriv2") );
            funcderiv3 = element.(strcat("Func",ab,"deriv3") );
            funcderiv4 = element.(strcat("Func",ab,"deriv4") );
            functdelay = element.(strcat("Func",ab,"TimeDelay") );
            turnidx = mod(turn, nsamples);

            t = t - functdelay;
            t2 = t * t;
            ampout = ampoutaux * ( ...
                func(turnidx) ...
                + funcderiv1(turnidx) * t ...
                + 0.5 * funcderiv2(turnidx) * t2 ...
                + 1.0 / 6.0 * funcderiv3(turnidx) * t2 * t ...
                + 1.0 / 24.0 * funcderiv4(turnidx) * t2 * t2 ...
            );
        else
            ampout = 0.0;
        end
      otherwise
        ampout = 0.0;
    end
end
