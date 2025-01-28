function elem=atvariablethinmultipole(fname, mode, varargin)
% ATVARIABLETHINMULTIPOLE Creates a variable thin multipole element
%
%  ATVARIABLETHINMULTIPOLE(FAMNAME,MODE)
%  ATVARIABLETHINMULTIPOLE(FAMNAME,MODE,PASSMETHOD,[KEY,VALUE]...)
%
%  INPUTS
%    FNAME          Family name 
%    MODE           Excitation mode: 'SINE', 'WHITENOISE' or 'ARBITRARY'.
%
%  OPTIONS
%    PASSMETHOD     Tracking function. Default: 'VariableThinMPolePass'
%
%  MORE OPTIONS (order does not matter)
%    AmplitudeA     Vector or scalar to define the excitation amplitude for
%                   PolynomA
%    AmplitudeB     Vector or scalar to define the excitation amplitude for
%                   PolynomB
%    FrequencyA     Frequency of SINE excitation for PolynomA. Unit Hz
%    FrequencyB     Frequency of SINE excitation for PolynomB. Unit Hz.
%    PhaseA         Phase of SINE excitation for PolynomA
%    PhaseB         Phase of SINE excitation for PolynomB
%    SinAabove      Limit the sin function to be above. Default -1.
%    SinBabove      Limit the sin function to be above. Default -1.
%    FuncA          ARBITRARY excitation turn-by-turn (tbt) list for
%                   PolynomA
%    FuncB          ARBITRARY excitation turn-by-turn (tbt) list for
%                   PolynomB
%    FuncAderiv1    ARBITRARY excitation tbt kick list for PolynomA 1st
%                   derivative wrt tau. Default: zeros(1,length(FUNC))
%    FuncBderiv1    ARBITRARY excitation tbt kick list for PolynomB 1st
%                   derivative wrt tau. Default: zeros(1,length(FUNC))
%    FuncAderiv2    ARBITRARY excitation tbt kick list for PolynomA 2nd
%                   derivative wrt tau. Default: zeros(1,length(FUNC))
%    FuncBderiv2    ARBITRARY excitation tbt kick list for PolynomB 2nd
%                   derivative wrt tau. Default: zeros(1,length(FUNC))
%    FuncAderiv3    ARBITRARY excitation tbt kick list for PolynomA 3rd
%                   derivative wrt tau. Default: zeros(1,length(FUNC))
%    FuncBderiv3    ARBITRARY excitation tbt kick list for PolynomB 3rd
%                   derivative wrt tau. Default: zeros(1,length(FUNC))
%    FuncAderiv4    ARBITRARY excitation tbt kick list for PolynomA 3rd
%                   derivative wrt tau. Default: zeros(1,length(FUNC))
%    FuncBderiv4    ARBITRARY excitation tbt kick list for PolynomB 3rd
%                   derivative wrt tau. Default: zeros(1,length(FUNC))
%    FuncATimeDelay TimeDelay to generate a small time offset on the
%                   function FUNC. It only has an effect if any of the
%                   derivatives is not zero. Default: 0.
%    FuncBTimeDelay TimeDelay to generate a small time offset on the
%                   function FUNC. It only has an effect if any of the
%                   derivatives is not zero. Default 0.
%    Periodic       If true the user input list is repeated. Default: false.
%    Ramps          Vector (t0, t1, t2, t3) in turn number to define the
%                   ramping of the excitation:
%                   * t<=t0: excitation amplitude is zero.
%                   * t0<t<=t1: excitation amplitude is linearly ramped up.
%                   * t1<t<=t2: excitation amplitude is constant.
%                   * t2<t<=t3: excitation amplitude is linearly ramped down.
%                   * t3<t:    excitation amplitude is zero.
%
%  OUTPUTS
%  1. ELEM - Structure with the AT element
%
%  NOTES
%    1. For all excitation modes at least one amplitude (A or B) is
%    required.
%    2. For SINE excitation modes the FREQUENCY is required.
%    3. For ARBITRARY excitation modes the FUNC is required.
%
%  EXAMPLES
%
% % Create a sinusoidal skew quadrupole with amplitude 0.0001  and frequency 1 kHz
% >> atvariablethinmultipole('ACSKEW','SINE','AmplitudeA',[0,1.e-4],'FrequencyA',1.e3);
%
% % Create a white noise horizontal dipole excitation of amplitude 0.1 mrad
% >> atvariablethinmultipole('ACKICK','WHITENOISE','AmplitudeB',1.e-4);
%
% % Create a vertical kick in the first turn and the opposite kick in the second
% % turn.
% >> funca = [1 -1 0];
% >> atvariablethinmultipole('CUSTOMFUNC','ARBITRARY','AmplitudeA',1e-4, ...
%                     'FuncA',funca,'Periodic',false);
%
%
% MORE DETAILS
%
% This function creates a thin multipole of any order (dipole kick, quadrupole,
% sextupole, etc.) and type (Normal or Skew) defined by the AmplitudeA and/or
% AmplitudeB components; the polynoms PolynomA and PolynomB are calculated on
% every turn depending on the chosen mode, and for some modes also on the
% particle time delay. All modes could be ramped.
%
% Keep in mind that as this element varies on every turn, and at the end of
% the tracking PolynomA and PolynomB are set to zero.
%
% Passing arrays of zeros as amplitude will initialize the MaxOrder to
% zero, and the polynom to a single zero.
%
% There are three different modes that could be set:
%   SINE (0), WHITENOISE (1) and ARBITRARY (2).
% when creating the element use 'SINE', 'WHITENOISE', or 'ARBITRARY'.
%
% The SINE mode requires amplitude and frequency for A and/or B.
% The jth component of the polynom (A or B) at the nth turn is given by:
%   amplitude(j)*sin[ TWOPI*frequency*(n*T0 + \tau_p) + phase],
% where T0 is the revolution period of the ideal ring, and \tau_p is the delay
% of the pth particle i.e. the sixth coordinate over the speed of light. Also,
% note that the position of the element on the ring has no effect, the phase
% could be used to add any delay due to the position along s.
% The following is an example of the SINE mode of an skew quad:
%     atvariablethinmultipole('VAR_SKEW','SINE',
%         'AmplitudeA',[0,skewa2],'FrequencyA',freqA,'PhaseA',phaseA)
% The values of the sin function could be limited to be above a defined
% threshold using ``Sin[AB]above``. For example, you could create a half-sin
% by setting ``Sin[AB]above`` to zero.
%
% The WHITENOISE mode requires the amplitude of either A or B. For example
%     atvariablethinmultipole('MYNOISE','WHITENOISE',
%         'AmplitudeA',[noiseA1])
% creates a gaussian vertical noise of amplitude noiseA1. The gaussian
% distribution is generated with zero-mean and one standard deviation from
% a pseudo-random stream pcg32. The pcg32 seed is fixed by the tracking
% function, therefore using the same stream on all trackings (sequencial or
% parallel). See https://github.com/atcollab/at/discussions/879 for more
% details on the pseudo random stream.
%
% The ARBITRARY mode requires the definition of a custom discrete function.
% The user defines the value of the function and its Taylor expansion with
% respect to \tau up to fourth order.
%     value = f(turn) + f'(turn)*tau + 0.5*f''(turn)*tau**2
%             + 1/6*f'''(turn)*tau**3 + 1/24*f''''(turn)*tau**4
% where f is an array of values, f',f'',f''',f'''', are arrays containing
% the derivatives wrt \tau, and \tau is the time delay of the particle, i.e.
% the the sixth coordinate divided by the speed of light.
% tau could be offset using FuncATimeDelay or FuncBTimeDelay.
%   tau <- tau - Func[AB]TimeDelay
% The function value is then multiplied by Amplitude A and/or B.
% For example, the following is a positive vertical kick in the first turn,
% negative on the second turn, and zero on the third turn.
%     atvariablethinmultipole('CUSTOMFUNC','ARBITRARY', ...
%         'AmplitudeA',1e-4,'FuncA',[1 -1 0],'Periodic',True);
% by default the array is assumed non periodic. The function has
% no effect on the particle in turns exceeding the function definition.
% If Periodic is set to True, the sequence is repeated.

%

% Input parser for option
[method,rsrc]=getargs(varargin,'VariableThinMPolePass', ...
    'check',@(arg) (ischar(arg) || isstring(arg)) && endsWith(arg,'Pass'));
[method,rsrc]                     = getoption(rsrc,'PassMethod',method);
[cl,rsrc]                         = getoption(rsrc,'Class','VariableMultipole');
rsrc                              = struct(rsrc{:});

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


    function rsrc = setsine(rsrc, ab)
        if ~isfield(rsrc,strcat('Frequency',ab))
            error(strcat('Please provide a value for Frequency',ab))
        end
        funcarg=strcat('Phase',ab);
        if ~isfield(rsrc,funcarg)
            rsrc.(funcarg) = 0;
        end
        funcarg=strcat('Sin',ab,'above');
        if ~isfield(rsrc,funcarg)
            rsrc.(funcarg) = -1;
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
            if strcmpi(mode,'SINE')
                rsrc = setsine(rsrc,ab);
            end
            if strcmpi(mode,'ARBITRARY')
                rsrc = setarb(rsrc,ab);
                if ~isfield(rsrc,'Periodic')
                    rsrc.Periodic = false;
                end
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
        mxab=max([mxa,mxb,1]);
        rsrc.MaxOrder=mxab-1;
        if isfield(rsrc,'AmplitudeA')
            rsrc.AmplitudeA(mxa+1:mxab)=0;
        end
        if isfield(rsrc,'AmplitudeB')
            rsrc.AmplitudeB(mxb+1:mxab)=0;
        end               
    end
    
end
