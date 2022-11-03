function elem=atvariablemultipole(fname,varargin)
%ATVARIABLEMULTIPOLE Creates a variable thin multipole element
%
%  ATVARIABLEMULTIPOLE(FAMNAME,MODE,PASSMETHOD,[KEY,VALUE]...)
%	
%  INPUTS
%    FNAME          Family name 
%    MODE           Excitation mode: 'SINE', 'WHITENOISE' or 'ARBITRARY'.
%                   Default: 'SINE'
%    PASSMETHOD     Tracking function. Default: 'VariableThinMPolePass'
%
%  OPTIONS (order does not matter)
%    AMPLITUDEA     Vector or scalar to define the excitation amplitude for
%                   PolynomA
%    AMPLITUDEB     Vector or scalar to define the excitation amplitude for
%                   PolynomA
%    FREQUENCYA     Frequency of SINE excitation for PolynomA
%    FREQUENCYB     Frequency of SINE excitation for PolynomB
%    PHASEA         Phase of SINE excitation for PolynomA
%    PHASEB         Phase of SINE excitation for PolynomB
%	 MAXORDER       Order of the multipole for a scalar amplitude
%    SEED           Input seed for the random number generator
%    FUNCA          ARBITRARY excitation turn-by-turn kick list for PolynomA
%    FUNCB          ARBITRARY excitation turn-by-turn kick list for PolynomB
%    PERIODIC       If true (default) the user input kick list is repeated
%    RAMPS          Vector (t0, t1, t2, t3) in turn number to define the ramping of the excitation
%                   * t<t0: excitation amlpitude is zero
%                   * t0<t<t1: exciation amplitude is linearly ramped up
%                   * t1<t<t2: exciation amplitude is constant             
%                   * t2<t<t3: exciation amplitude is linearly ramped down
%                   * t3<t: exciation amplitude is zero 
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
[mode,rsrc]=getargs(varargin,'SINE','check',@(arg) any(strcmpi(arg,{'SINE','WHITENOISE','ARBITRARY'})));
[method,rsrc]=getargs(rsrc,'VariableThinMPolePass','check',@(arg) (ischar(arg) || isstring(arg)) && endsWith(arg,'Pass'));
[mode,rsrc]                       = getoption(rsrc,'Mode',mode);
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
        rsrc.(strcat('NSamples',ab))=length(rsrc.(funcarg));
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
        else
            mxa=0;
        end
        if isfield(rsrc,'AmplitudeB')
            mxb=find(abs(rsrc.AmplitudeB)>0,1,'last');
        else
            mxb=0;
        end
        mxab=max([mxa,mxb,rsrc.MaxOrder-1]);
        rsrc.MaxOrder=mxab-1;
        if isfield(rsrc,'AmplitudeA')
            rsrc.AmplitudeA(mxa+1:mxab)=0;
        end
        if isfield(rsrc,'AmplitudeB')
            rsrc.AmplitudeB(mxb+1:mxab)=0;
        end               
    end
    
end
