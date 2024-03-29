function newring=atfittune(ring,varargin)
%ATFITTUNE Fit linear tunes by scaling 2 quadrupole families
%
% NEWRING = ATFITTUNE(RING,NEWTUNES,QUADFAMILY1,QUADFAMILY2)
% NEWRING = ATFITTUNE(RING,DPP,NEWTUNES,QUADFAMILY1,QUADFAMILY2)
%
%RING:          Cell array
%DPP:           Optional momentum deviation (default 0)
%NEWTUNES:      Desired tune values
%QUADFAMILY1:   1st quadrupole family
%QUADFAMILY2:   2nd quadrupole family
%
%QUADFAMILY may be:
%   string: Family name
%   logical array: mask of selected elements in RING
%   Numeric array: list of selected elements in RING
%   Cell array: All elements selected by each cell
%
%NEWRING = ATFITTUNE(RING,...,'UseIntegerPart') With this flag, the
%   function fits the tunes to the total values of NEWTUNES, including
%   the integer part.
%   With this option the function is substantially slower!
%
%NEWRING = ATFITTUNE(RING,...,'KStep',kstep)
%   kstep is the quadrupole strength applied to build the jacobian [m^-2].
%   Default: 1.0e-6
%
%NEWRING = ATFITTUNE(RING,...,'dp',DP)
%   Specify off-momentum if radiation is OFF (default 0)
%
%NEWRING = ATFITTUNE(RING,...,'dct',DCT)
%   Specify path lengthening if radiation is OFF (default 0)
%
% See also ATFITCHROM

newargs=getdparg(varargin);
newring=wrapper6d(ring,@fittune,newargs{:});

    function newring=fittune(ring,~,varargin)
        [UseIntegerPart,varargs]=getflag(varargin,'UseIntegerPart');
        [delta,varargs]=getoption(varargs,'KStep',1.0e-6);
        [newtunes,famname1,famname2,varargs]=getargs(varargs,[],[],[]);

        idx1=varelem(ring,famname1);
        idx2=varelem(ring,famname2);

        kl1=atgetfieldvalues(ring(idx1),'PolynomB',{2});
        kl2=atgetfieldvalues(ring(idx2),'PolynomB',{2});

        if UseIntegerPart
            allpos=1:length(ring)+1;
            gettune = @getinttune;
        else
            if any(newtunes>=1)
                warning('AT:FitTune','The integer part of the tunes is ignored unless you use the ''UseIntegerPart'' flag');
                newtunes=newtunes-floor(newtunes);
            end
            gettune = @getfractune;
        end

        % Compute initial tunes before fitting
        tunes = gettune(ring,varargs{:});

        % Take Derivative
        tunes1 = gettune(setqp(ring,idx1,kl1,delta),varargs{:});
        tunes2 = gettune(setqp(ring,idx2,kl2,delta),varargs{:});

        %Construct the Jacobian
        J = [tunes1-tunes tunes2-tunes]/delta;
        dK = J\(newtunes(:)-tunes);

        % Apply new strengths
        newring = setqp(ring,idx1,kl1,dK(1));
        newring = setqp(newring,idx2,kl2,dK(2));

        function ring2=setqp(ring,idx,k0,delta)
            k=k0*(1+delta);
            ring2=ring;
            ring2(idx)=atsetfieldvalues(atsetfieldvalues(ring2(idx),'K',k),'PolynomB',{2},k);
        end

        function res=varelem(ring,arg)
            if islogical(arg)
                res=arg;
            elseif isnumeric(arg)
                res=false(size(ring));
                res(arg)=true;
            elseif ischar(arg)
                res=atgetcells(ring,'FamName',arg);
            elseif iscell(arg)
                res=false(size(ring));
                for i=1:length(arg)
                    res=res|varelem(ring,arg{i});
                end
            else
                error('AT:GetElemList:WrongArg','Cannot parse argument');
            end
        end

        function tun=getfractune(ring,varargin)
            tun3=tunechrom(ring,varargin{:});
            tun=tun3(1:2)';
        end
        
        function tun=getinttune(ring,varargin)
            [ringdata,elemdata] = atlinopt6(ring,allpos,varargin{:}); %#ok<ASGLU>
            tun=elemdata(end).mu(1:2)'/2/pi;
        end
    end
end
