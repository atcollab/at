function elem=atvariablemultipole(fname, mode, varargin)
% ATVARIABLEMULTIPOLE Creates a variable thin multipole element
%
%  ATVARIABLEMULTIPOLE(FAMNAME,MODE)
%  ATVARIABLEMULTIPOLE(FAMNAME,MODE,PASSMETHOD,[KEY,VALUE]...)
%
% This function creates a thin multipole of order and type defined by the
% amplitude components; the polynoms PolynomA and PolynomB are calculated
% on every turn depending on the chosen mode.
%
% Keep in mind that this element varies on every turn, therefore, any ring
% containing a variable element may change after tracking n turns.
%
% There are three different modes implemented: SINE, WHITENOISE and ARBITRARY.
%
% The SINE mode requires amplitude, frequency and phase of at least one of the
% two polynoms A or B. The j-th component of the polynom on the n-th turn
% is given by:
%   amplitude_j*sin[ 2\pi*frequency*(nth_turn*T0 + c\tau_k) + phase],
% where T0 is the revolution period of the ideal ring, and c\tau_k is the
% delay of the kth particle i.e. the sixth coordinate over the speed of light.
% Also, note that the position of the element on the ring has no effect,
% the phase should be used to add any delay due to the position along s.
% The following is an example of the SINE mode of an skew quad:
%
%     varskew = ATVARIABLEMULTIPOLE('VAR_SKEW','SINE', ...
%         AmplitudeA=[0,skewa2],FrequencyA=freqA,PhaseA=phaseA)
%
% The WHITENOISE mode requires the amplitude.
% THe ARBITRARY mode requires the amplitude
%
%
%  INPUTS
%    FNAME          Family name 
%    MODE           Excitation mode: 'SINE', 'WHITENOISE' or 'ARBITRARY'.
%
%  OPTIONS (order does not matter)
%    PASSMETHOD     Tracking function. Default: 'VariableThinMPolePass'
%    AMPLITUDEA     Vector or scalar to define the excitation amplitude for
%                   PolynomA
%    AMPLITUDEB     Vector or scalar to define the excitation amplitude for
%                   PolynomB
%    FREQUENCYA     Frequency of SINE excitation for PolynomA
%    FREQUENCYB     Frequency of SINE excitation for PolynomB
%    PHASEA         Phase of SINE excitation for PolynomA
%    PHASEB         Phase of SINE excitation for PolynomB
%    SEED           Input seed for the random number generator
%    FUNCA          ARBITRARY excitation turn-by-turn (tbt) kick list for PolynomA
%    FUNCB          ARBITRARY excitation turn-by-turn (tbt) kick list for PolynomB
%    FUNCADERIV1    ARBITRARY excitation tbt kick list for PolynomA 1st
%                   derivative wrt tau
%    FUNCBDERIV1    ARBITRARY excitation tbt kick list for PolynomB 1st
%                   derivative wrt tau
%    FUNCADERIV2    ARBITRARY excitation tbt kick list for PolynomA 2nd
%                   derivative wrt tau
%    FUNCBDERIV2    ARBITRARY excitation tbt kick list for PolynomB 2nd
%                   derivative wrt tau
%    FUNCADERIV3    ARBITRARY excitation tbt kick list for PolynomA 3rd
%                   derivative wrt tau
%    FUNCBDERIV3    ARBITRARY excitation tbt kick list for PolynomB 3rd
%                   derivative wrt tau
%    FUNCTIMEDELAY  TimeDelay to generate a small time offset on the
%                   function FUNC. It only has an effect if any of the
%                   derivatives is not zero.
%    PERIODIC       If true (default) the user input kick list is repeated
%    RAMPS          Vector (t0, t1, t2, t3) in turn number to define the ramping of the excitation
%                   * t<t0: excitation amplitude is zero
%                   * t0<t<t1: excitation amplitude is linearly ramped up
%                   * t1<t<t2: excitation amplitude is constant
%                   * t2<t<t3: excitation amplitude is linearly ramped down
%                   * t3<t:    excitation amplitude is zero
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
%
% % Create a sinusoidal dipole with amplitude 0.1 mrad and frequency 1 kHz
% >> atvariablemultipole('ACM','SINE','AmplitudeB',1.e-4,'FrequencyB',1.e3);
%
% % Create a white noise dipole excitation of amplitude 0.1 mrad
% >> atvariablemultipole('ACM','WHITENOISE','AmplitudeB',1.e-4);

% Input parser for option
[method,rsrc]=getargs(varargin,'VariableThinMPolePass', ...
    'check',@(arg) (ischar(arg) || isstring(arg)) && endsWith(arg,'Pass'));
[method,rsrc]                     = getoption(rsrc,'PassMethod',method);
[cl,rsrc]                         = getoption(rsrc,'Class','VariableMultipole');
[maxorder,rsrc]                   = getoption(rsrc,'MaxOrder',0);
rsrc                              = struct(rsrc{:});
rsrc.MaxOrder                     = maxorder;

m=struct('SINE',0,'WHITENOISE',1,'ARBITRARY',2);

if ~any(isfield(rsrc,{'AmplitudeA','AmplitudeB'}))
    error("Please provide at least one amplitude for A or B")
end
rsrc = setparams(rsrc,mode,'A');
rsrc = setparams(rsrc,mode,'B');
rsrc = setmaxorder(rsrc);

% Build the element
% rsrc =namedargs2cell(rsrc);   % introduced in R2019b
rsrc=reshape([fieldnames(rsrc) struct2cell(rsrc)]',1,[]);
elem=atbaselem(fname,method,'Class',cl,'Length',0,'Mode',m.(upper(mode)),...
               'PolynomA',[],'PolynomB',[],rsrc{:});


    function setsine(rsrc, ab)
        if ~isfield(rsrc,strcat('Frequency',ab))
            error(strcat('Please provide a value for Frequency',ab))
        end
    end

    function rsrc = setarb(rsrc, ab)
        funcarg=strcat('Func',ab);
        if ~isfield(rsrc,funcarg)
            error(strcat('Please provide a value for Func',ab))
        end
        nsamples = length(rsrc.(funcarg));
        rsrc.(strcat('NSamples',ab)) = nsamples;
        for i = 1:4
            funcarg=strcat('Func',ab,'deriv',num2str(i));
            if ~isfield(rsrc,funcarg)
                rsrc.(funcarg) = zeros(1,nsamples);
            end
        end
        funcarg = strcat('Func',ab,'TimeDelay');
        if ~isfield(rsrc,funcarg)
            rsrc.(funcarg) = 0;
        end
    end

    function rsrc = setparams(rsrc,mode,ab)
        amplarg=strcat('Amplitude',ab);
        if isfield(rsrc,amplarg)
            amp=rsrc.(amplarg);
            if isscalar(amp)
                rsrc.(amplarg)=[zeros(1,rsrc.MaxOrder) amp];
            end
            if strcmpi(mode,'SINE')
                setsine(rsrc,ab);
            end
            if strcmpi(mode,'ARBITRARY')
                rsrc = setarb(rsrc,ab);
            end        
        end
    end

    function rsrc = setmaxorder(rsrc)
        if isfield(rsrc,'AmplitudeA')
            mxa=find(abs(rsrc.AmplitudeA)>0,1,'last');
            if isempty(mxa)
                mxa=1;
            end
        else
            mxa=0;
        end
        if isfield(rsrc,'AmplitudeB')
            mxb=find(abs(rsrc.AmplitudeB)>0,1,'last');
            if isempty(mxb)
                mxb=1;
            end
        else
            mxb=0;
        end
        mxab=max([mxa,mxb,rsrc.MaxOrder+1]);
        rsrc.MaxOrder=mxab-1;
        if isfield(rsrc,'AmplitudeA')
            rsrc.AmplitudeA(mxa+1:mxab)=0;
        end
        if isfield(rsrc,'AmplitudeB')
            rsrc.AmplitudeB(mxb+1:mxab)=0;
        end               
    end
    
end
