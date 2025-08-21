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
%    BufferSizeA    Set the buffer length in WHITENOSE mode.
%    BufferSizeB    Set the buffer length in WHITENOSE mode.
%    FuncA          ARBITRARY excitation turn-by-turn (tbt) list for
%                   PolynomA
%    FuncB          ARBITRARY excitation turn-by-turn (tbt) list for
%                   PolynomB
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
% This Class creates a thin multipole of any ord
% sextupole, etc.) and type (Normal or Skew) def
% AmplitudeB components; the polynoms PolynomA a
% on every turn depending on the chosen mode, an
% particle time delay. All modes could be ramped
% Default pass method: ``VariableThinMPolePass``
%
% Keep in mind that as this element varies on ev
% the tracking PolynomA and PolynomB are set to 
% with optics calculations. Therefore, it is pre
% to IdentityPass when the element should not be
%
% There are three different modes that could be 
% :py:attr:`.ACMode.SINE`: sine function
% :py:attr:`.ACMode.WHITENOISE`: gaussian white 
% :py:attr:`.ACMode.ARBITRARY`: user defined tur
% SINE = 0, WHITENOISE = 1 and ARBITRARY = 2. Se
%
% For example, use at.ACMode.SINE or mode = 0 to
%
% The **SINE** mode requires amplitude and frequ
% The value of the jth component of the polynom 
% is given by
%
% Amplitude[j]*sin[TWOPI*frequency*(n*T0 + \tau_
%
% where T0 is the revolution period of the ideal
% of the pth particle i.e. the sixth coordinate 
% note that the position of the element on the r
% could be used to add any delay due to the posi
% an example of the SINE mode of an skew quad
%
% eleskew = at.VariableThinMultipole('VAR_SKEW',
% AmplitudeA=[0,skewa2],FrequencyA=freqA,PhaseA=
%
% The values of the sin function could be limite
% threshold using ``Sin[AB]above``. For example,
% half-sin by setting ``Sin[AB]above`` to zero. 
% negative half-sin function by setting the ampl
%
% The **WHITENOISE** mode requires the amplitude
% distribution is generated with zero-mean and o
% a pseudo-random stream pcg32. The pcg32 seed i
% function, therefore using the same stream on a
% parallel). See https://github.com/atcollab/at/
% details on the pseudo random stream. For examp
%
% elenoise = at.VariableThinMultipole('MYNOISE',
% AmplitudeA=[noiseA1])
%
% creates a vertical kick as gaussian noise of a
%
% Additionally, one could use the parameter Buff
% to create buffer to store the random numbers u
% during tracking. The buffer could be smaller o
% If the buffer is smaller, it will store from s
% completely the buffer. If the buffer is larger
% to zero.
%
% Care should be taken when changing the buffer 
% done is at the moment of initialization, for e
% from a file. Buffer and buffer size could be r
% ensure that the actual buffer length and buffe
%
% The **ARBITRARY** mode requires the definition
% to be sampled at every turn. The function and 
% respect to \tau up to any given order is
%
% value = f(turn) + f'(turn)*tau + 0.5*f''(turn)
% + 1/6*f'''(turn)*tau**3 + 1/24*f''''(turn)*tau
%
% f is an array of values, f',f'',f''',f'''', ar
% the derivatives wrt \tau, and \tau is the time
% the the sixth coordinate divided by the speed 
% function func is a 2D-array with columns corre
% derivatives, and rows to turns. For example, a
% f(1)=1 with positive derivative f'(1)=0.1 foll
% second turn f(2)=-1 with negative derivative f
% Func=[[1,  -1];
%       [0.1,-0.2]].
% The time tau could be offset using ``FuncATime
%
% tau = tau - Func[AB]TimeDelay
%
% The function value is then **multiplied by Amp
% Use FuncA or FuncB, and AmplitudeA or Amplitud
% For example, the following is a positive verti
% negative on the second turn, and zero on the t
%
% elesinglekick = at.VariableThinMultipole('CUST
% AmplitudeA=1e-4,FuncA=[1,-1,0],Periodic=True)
%
% As already mentioned, one could include the de
% from a Taylor expansion in a 2D array. First r
% every turn, second row for the values of f'' o
%
% By default the array is assumed non periodic, 
% on the particle in turns exceeding the functio
% ``Periodic`` is set to True, the sequence is r
%
% One could use the method inspect_polynom_value
% used in every turn.
%
% Any mode could be ramped. The ramp is defined 
%     * ``t<=t0``: excitation amplitude is zero.
%     * ``t0<t<=t1``: excitation amplitude is li
%     * ``t1<t<=t2``: excitation amplitude is co
%     * ``t2<t<=t3``: excitation amplitude is li
%     * ``t3<t``: excitation amplitude is zero.
% See also:
%
%  ATINSPECTVARIABLETHINMULTIPOLE

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
        [ktaylor, nsamples] = size(rsrc.(funcarg));
        rsrc.(strcat('NSamples',ab)) = nsamples;
        rsrc.(strcat('Ktaylor', ab)) = ktaylor;
        funcarg = strcat('Func',ab,'TimeDelay');
        if ~isfield(rsrc,funcarg)
            rsrc.(funcarg) = 0;
        end
    end

    function rsrc = setwhitenoise(rsrc, ab)
      funcarg = strcat('BufferSize',ab);
      if ~isfield(rsrc,funcarg)
          rsrc.(funcarg) = 0;
      end
      funcarg = strcat('Buffer',ab);
      if ~isfield(rsrc,funcarg)
          bfsize = rsrc.(strcat('BufferSize',ab));
          rsrc.(funcarg) = zeros(1,bfsize);
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
            if strcmpi(mode,'WHITENOISE')
              rsrc = setwhitenoise(rsrc,ab);
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
