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
%    MAXORDER       Order of the multipole for a scalar amplitude
%    SEED           Input seed for the random number generator
%    FUNCA          ARBITRARY excitation turn-by-turn kick array for PolynomA
%    FUNCB          ARBITRARY excitation turn-by-turn kick array for PolynomB
%    FUNCDELAY      Value to substract from the particles 6th coordinate
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
%       required.
%    2. For ARBITRARY excitation modes the FUNC corresponding to the input
%       AMPLITUDE is required
%    3. In ARBITRARY mode the seed is fixed by the tracking function, and
%       it is common to all threads. See ringpass, linepass and elempass.
%    4. Func(A,B) could be an array of size (m,n) with n coefficients in the first
%       row for the function over n turns, and other m-1 rows with higher order
%       derivatives with respect to ctau
%       i.e. on the kth turn the ith component of the Polynom(A/B) seen by a
%       particle with 6th coordinate ctau is calculated as,
%          amp(A/B)[i] * func(0,k) + (ctau-delay) * func(1,k) + ...
%                                      + (ctau-delay)**m * func(m,k)
%
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
% % Create a pulse that decreses to 10% on the second turn and dissappears afterwards.
% >> atvariablethimultipole('PULSE','ARBITRARY','AmplitudeB',1e-3, ...
%        'FuncB',[1 0.1])
%
% % Delay the pulse to match the particle delay ctau=0.2 [m], and include the first order
%   derivative with respect to ctau
% >> atvariablethimultipole('PULSE','ARBITRARY','AmplitudeB',1e-3, ...
%        'FuncB',[1 0.1;der1turn0 der1turn1],'FuncDelay',0.2)


% Input parser for option

[modename, rsrc] = getargs(varargin,'SINE', ...
                   'check',@(arg) any(strcmpi(arg,{'SINE','WHITENOISE','ARBITRARY'})));
[modename, rsrc] = getoption(rsrc,'ModeName',modename);
[~, rsrc] = getoption(rsrc,'Mode',2); % remove Mode, it is set by ModeName
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

if modename == "ARBITRARY"
  if isfield(rsrc,'AmplitudeA') && ~isfield(rsrc,'FuncA')
    error("This mode requires matching Func(A/B) and Amplitude(A/B)");
  end
  if isfield(rsrc,'AmplitudeB') && ~isfield(rsrc,'FuncB')
    error("This mode requires matching Func(A/B) and Amplitude(A/B)");
  end
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

    function rsrc = setarb(rsrc, ab)
        funcarg=strcat('Func',ab);
        if isfield(rsrc,funcarg)
          [rows, nsamples] = size(rsrc.(funcarg));
          rsrc.(strcat('NSamples',ab)) = nsamples;
          rsrc.(strcat('Dorder', ab)) = rows-1;
        end
    end

    function rsrc = setparams(rsrc,modename,ab)
        amplarg=strcat('Amplitude',ab);
        if isfield(rsrc,amplarg)
            switch modename
                case "SINE"
                    rsrc = setsine(rsrc,ab);
                case "ARBITRARY"
                    rsrc = setarb(rsrc,ab);
                    if ~isfield(rsrc,'Periodic')
                       rsrc.Periodic = 0;
                    end
                    if ~isfield(rsrc,'FuncDelay')
                       rsrc.FuncDelay = 0;
                    end
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
