function elem=atvariablethinmultipole(fname,varargin)
%ATVARIABLETHINMULTIPOLE Creates a variable thin multipole element
%
%  ATVARIABLETHINMULTIPOLE(FAMNAME,[KEY,VALUE]...)
%  ATVARIABLETHINMULTIPOLE(FAMNAME,MODENAME,[KEY,VALUE]...)
%  ATVARIABLETHINMULTIPOLE(FAMNAME,MODENAME,PASSMETHOD,[KEY,VALUE]...)
%
%  INPUTS
%    FNAME          Family name
%
%  OPTIONS (order does not matter)
%    MODENAME       'SINE', 'WHITENOISE' or 'ARBITRARY'.
%                   Default: 'SINE'
%    PASSMETHOD     Tracking function. Default: 'VariableThinMPolePass'
%    AMPLITUDEA     Vector or scalar to define the excitation amplitude for
%                   PolynomA
%    AMPLITUDEB     Vector or scalar to define the excitation amplitude for
%                   PolynomA
%    FREQUENCYA     Frequency of SINE excitation for PolynomA
%    FREQUENCYB     Frequency of SINE excitation for PolynomB
%    PHASEA         Phase of SINE excitation for PolynomA
%    PHASEB         Phase of SINE excitation for PolynomB
%    SINMIN         Sine function min limit. Default -1.1
%    SINMAX         Sine function max limit. Default +1.1
%    BUFFERSIZE     Set the buffer length in WHITENOISE mode.
%    BUFFEROFFSET   Offset in turns for the first element of the buffer
%                         in WHITENOISE mode.
%    MAXORDER       Order of the multipole for a scalar amplitude
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
%    required.
%    2. For SINE excitation modes the FREQUENCY corresponding to the input
%    AMPLITUDE is required
%    3. For ARBITRARY excitation modes the FUNC corresponding to the input
%    AMPLITUDE is required
%    4. In ARBITRARY mode the seed is fixed by the tracking function, and
%    it is common to all threads. See at.track.
%
%  EXAMPLES
%
% % Create a sinusoidal dipole with amplitude 0.1 mrad and frequency 1 kHz
% >> atvariablethinmultipole('ACM','SINE','AmplitudeB',1.e-4,'FrequencyB',1.e3);
%
% % Create a white noise dipole excitation of amplitude 0.1 mrad
% >> atvariablethinmultipole('ACM','WHITENOISE','AmplitudeB',1.e-4);
%
% % Create a half-sine function
% >> atvariablethimultipole('HSINE','SINE','AmplitudeB',1e-3,'FrequencyB',100,'Sinmin',0)
%
% % Create a sine saturation
% >> atvariablethimultipole('HSINE','SINE','AmplitudeB',1e-3,'FrequencyB',100,'Sinmax',0.9)
%
% % Create whitenoise and store the random values starting from third turn for 10 turns.
% >> atvariablethinmultipole('v',"WHITENOISE",'AmplitudeB',1e-3,'BufferSize',10,'BufferOffset',2);

% Input parser for option

[modename, rsrc] = getargs(varargin,'SINE', ...
                   'check',@(arg) any(strcmpi(arg,{'SINE','WHITENOISE','ARBITRARY'})));
[modename, rsrc] = getoption(rsrc,'ModeName',modename);
if ~any(strcmpi(modename,{'SINE','WHITENOISE','ARBITRARY'}))
  error("ModeName should be 'SINE', 'WHITENOISE' or 'ARBITRARY'");
end
[method,rsrc]   = getargs(rsrc,'VariableThinMPolePass', ...
                  'check',@(arg) (ischar(arg) || isstring(arg)) && endsWith(arg,'Pass'));
[method,rsrc]   = getoption(rsrc,'PassMethod',method);
[cl,rsrc]       = getoption(rsrc,'Class','VariableThinMultipole');
[maxorder,rsrc] = getoption(rsrc,'MaxOrder',0);
rsrc            = struct(rsrc{:});
rsrc.MaxOrder   = maxorder;

if ~any(isfield(rsrc,{'AmplitudeA','AmplitudeB'}))
    rsrc.AmplitudeB = 0;
end

rsrc = setparams(rsrc,modename,'A');
rsrc = setparams(rsrc,modename,'B');
rsrc = setmaxorder(rsrc);

m=struct('SINE',0,'WHITENOISE',1,'ARBITRARY',2);

% Build the element
% rsrc =namedargs2cell(rsrc);   % introduced in R2019b
rsrc=reshape([fieldnames(rsrc) struct2cell(rsrc)]',1,[]);
elem=atbaselem(fname,method,'Class',cl,'Length',0,'Mode',m.(modename),...
               'ModeName',modename,'PolynomA',[],'PolynomB',[],rsrc{:});


    function rsrc = setsine(rsrc, ab)
        funcarg = strcat('Frequency',ab);
        if ~isfield(rsrc,funcarg)
            rsrc.(funcarg) = 0;
        end
        funcarg=strcat('Phase',ab);
        if ~isfield(rsrc,funcarg)
            rsrc.(funcarg) = 0;
        end
        funcarg="Sinmin";
        if ~isfield(rsrc,funcarg)
            rsrc.(funcarg) = -1.1;
        end
        funcarg="Sinmax";
        if ~isfield(rsrc,funcarg)
            rsrc.(funcarg) = 1.1;
        end

    end

    function rsrc = setwhitenoise(rsrc, ab)
      funcarg = 'BufferSize';
      if ~isfield(rsrc,funcarg)
          rsrc.(funcarg) = 0;
      end
      funcarg = 'BufferOffset';
      if ~isfield(rsrc,funcarg)
          rsrc.(funcarg) = 0;
      end
      funcarg = strcat('Buffer',ab);
      if ~isfield(rsrc,funcarg)
          bfsize = rsrc.('BufferSize');
          rsrc.(funcarg) = zeros(1,bfsize);
      end
    end

    function rsrc = setarb(rsrc, ab)
        funcarg=strcat('Func',ab);
        if ~isfield(rsrc,funcarg)
            error(strcat('Please provide a value for Func',ab))
        end
        rsrc.(strcat('NSamples',ab))=length(rsrc.(funcarg));
    end

    function rsrc = setparams(rsrc,modename,ab)
        amplarg=strcat('Amplitude',ab);
        if isfield(rsrc,amplarg)
            switch modename
                case "SINE"
                    rsrc = setsine(rsrc,ab);
                case "ARBITRARY"
                    rsrc = setarb(rsrc,ab);
                case "WHITENOISE"
                    rsrc = setwhitenoise(rsrc,ab);
            end
        end
    end

    function rsrc = setmaxorder(rsrc)
        mxab = [0 0];
        thefields = {'AmplitudeA','AmplitudeB'};
        for i = 1:2
          ampab = thefields{i};
          tmp = 0;
          if isfield(rsrc,ampab)
              tmp=find(abs(rsrc.(ampab))>0,1,'last');
              if isempty(tmp), tmp=1; end
          end
          mxab(i) = tmp;
        end
        maxamp=max([mxab,rsrc.MaxOrder+1]);
        rsrc.MaxOrder=maxamp-1;
        for i = 1:2
          ampab = thefields{i};
          if isfield(rsrc,ampab)
            rsrc.(ampab)(mxab(i)+1:maxamp)=0;
          end
        end
    end
end
