function elem=atvariablemultipole(fname,varargin)
%ATVARIABLEMULTIPOLE Creates a variable thin multipole element
%
%  ATVARIABLEMULTIPOLE(FAMNAME,MODE,PASSMETHOD)
%	
%  INPUTS
%  1. FNAME      - Family name 
%  2. MODE       - Excitation mode: SINE, WHITENOISE, ARBITRARY
%  5. PASSMETHOD - Tracking function. Defaults to 'VariableThinMPolePass'
%
%  OPTIONS (order does not matter)
%    AMPLITUDEA  -	Vector or scalar to define the excitation amplitude for
%                   PolynomA
%    AMPLITUDEB  -	Vector or scalar to define the excitation amplitude for
%                   PolynomA
%    FREQUENCYA  -	Frequency of SINE excitation for PolynomA
%    FREQUENCYB  -	Frequency of SINE excitation for PolynomB
%    PHASEA      -	Phase of SINE excitation for PolynomA
%    PHASEB      -	Phase of SINE excitation for PolynomB
%	 MAXORDER    -  Order of the multipole for a scalar amplitude
%	 SEED        -  Seed of the random number generator for the white noise
%                   excitation
%    FUNCA      -	ARBITRAY excitation turn-by-turn kick list for PolynomA
%    FUNCB      -	ARBITRAY excitation turn-by-turn kick list for PolynomB
%    RAMPS      -	Vector of length 4 to determine the beginning and end
%                   of the  linear ramp up and down
%
%  OUTPUTS
%  1. ELEM - Structure with the AT element
%
%  NOTES
%    1. For all excitation modes at least one amplitude (A or B) is
%    required. The default excitation is SINE
%    2. For SINE excitation modes the FREQUENCY corresponding to the input
%    AMPLITUDE is required
%    3. For ARBITRARY excitation modes the FUNC corresponding to the input
%    AMPLITUDE is required
%
%  EXAMPLES
%    1. atvariablemultipole(famname, mode, passmethod,'fieldname1',value1,...)
%   each pair {'fieldname',value} is added to the element
%


% Input parser for option
[rsrc,mode,method] = decodeatargs({'SINE','VariableThinMPolePass'},varargin);
[mode,rsrc]                       = getoption(rsrc,'Mode',mode);
[method,rsrc]                     = getoption(rsrc,'PassMethod',method);
[cl,rsrc]                         = getoption(rsrc,'Class','VariableMultipole');
[maxorder,rsrc]                   = getoption(rsrc,'MaxOrder',0);
rsrc                              = struct(rsrc{:});
rsrc.MaxOrder                     = maxorder;

keys = {'SINE','WHITENOISE','ARBITRARY'};
vals = [0 1 2];
m = containers.Map(keys,vals);

a = isfield(rsrc,'AmplitudeA');
b = isfield(rsrc,'AmplitudeB');

if ~a && ~b
    error("Please provide at least one amplitude for A or B")
end
rsrc = setparams(rsrc,mode,'A');
rsrc = setparams(rsrc,mode,'B');
rsrc = setmaxorder(rsrc);

% Build the element
rsrc =namedargs2cell(rsrc);
elem=atbaselem(fname,method,'Class',cl,'Length',0,'Mode',m(mode),...
               'PolynomA',[],'PolynomB',[],rsrc{:});


    function setsine(rsrc, ab)
        if ~isfield(rsrc,{strcat('Frequency',ab)})
            error(strcat('Please provide a value for Frequency',ab))
        end
    end

    function rsrc = setarb(rsrc, ab)
        if ~isfield(rsrc,{strcat('Func',ab)})
            error(strcat('Please provide a value for Func',ab))
        end
        nsamp=length(rsrc.(strcat('Func',ab)));
        rsrc.(strcat('NSamples',ab))=nsamp;
    end

    function rsrc = setparams(rsrc,mode,ab)
        if isfield(rsrc,strcat('Amplitude',ab))
            amp=rsrc.(strcat('Amplitude',ab));
            if isscalar(amp)
                amp=[zeros(1,rsrc.MaxOrder) amp];
            end
            rsrc.(strcat('Amplitude',ab))=amp;
            if strcmp(mode,'SINE')
                setsine(rsrc,ab);
            end
            if strcmp(mode,'ARBITRARY')
                rsrc = setarb(rsrc,ab);
            end        
        end
    end

    function rsrc = setmaxorder(rsrc)
        mxa=0;
        mxb=0;
        if isfield(rsrc,'AmplitudeA')
            mxa=find(abs(rsrc.AmplitudeA)>0,1,'last')-1;
        end
        if isfield(rsrc,'AmplitudeB')
            mxb=find(abs(rsrc.AmplitudeB)>0,1,'last')-1;
        end
        rsrc.MaxOrder=max([mxa,mxb]);
        if isfield(rsrc,'AmplitudeA')
            delta=rsrc.MaxOrder-length(rsrc.AmplitudeA)+1;
            if delta>0
                rsrc.AmplitudeA=[rsrc.AmplitudeA zeros(1,delta)];
            end
        end
        if isfield(rsrc,'AmplitudeB')
            delta=rsrc.MaxOrder-length(rsrc.AmplitudeB)+1;
            if delta>0
                rsrc.AmplitudeB=[rsrc.AmplitudeB zeros(1,delta)];
            end
        end               
    end
    
end
